import os, logging, time
import numpy as np

from scipy.stats import truncnorm
from scipy.special import binom as sp_binom, erfinv
from scipy.special import gammaln#factorial as sp_factorial

from collections import Counter
import quadpy

from matplotlib import pyplot as plt

import ultranest
import ultranest.stepsampler
from ultranest.plot import cornerplot
from dynesty import DynamicNestedSampler, pool as dypool, utils, plotting as dyplot


logging.basicConfig(level=logging.INFO)



def p_nu_single(NU,gamma,delta,nu_max,log=False):

    NU_mask = NU > 0
    p = np.zeros_like(NU)
    NU_scaled = NU[NU_mask] / nu_max
    if log:
        return - np.log( nu_max / gamma * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) - delta**2 / 2 + \
        ( gamma**2 - 1 ) * np.log( NU / nu_max ) + \
        np.log( np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max ) ) ) )
    else:
        p[NU_mask] = gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU_scaled ) ) ) * \
            np.exp( - delta**2/2.) * ( NU_scaled )**(gamma**2 - 1) * \
            np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU_scaled ) ) )
        p[~NU_mask] = 0 if gamma > 1 else np.inf
        return p

def p_nu(NU,p,two_pop=False):

    if two_pop:
        return (p['weight_dark'] * p_nu_single(NU,p['gamma_dark'],p['delta_dark'],p['nu_max']) + \
        (1-p['weight_dark']) * p_nu_single(NU,p['gamma'],p['delta'],p['nu_max']))
    else:
        return p_nu_single(NU,p['gamma'],p['delta'],p['nu_max'])

def poisson_spikes(nu,N_AP,T,log=False):
    ## using the gamma-function to obtain log(N!) = gammaln(N+1)
    if log:
        return N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T
    else:
        return np.exp(N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T)


def f(NU,p,N_AP,T,zero=False,two_pop=False,log=False):
    '''
        calculates the probability to observe N_AP action potentials in any neuron, 
        given the underlying firing rate distribution and poisson firing
    '''

    if zero:
        return p_nu(NU,p,two_pop=two_pop) * np.exp(-NU*T)
    else:
        return p_nu(NU,p,two_pop=two_pop) * poisson_spikes(NU,N_AP[:,np.newaxis],T,log=log)



class BayesModel:


    def __init__(self,modelParams,mode='rates',logLevel=logging.ERROR):
        '''
            initialize the model by setting the priors and preparing the data
        '''

        self.log = logging.getLogger("nestLogger")
        self.set_logLevel(logLevel)

        os.environ['MKL_NUM_THREADS'] = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        os.environ['OMP_NUM_THREADS'] = '1'


        self.prepare_data(modelParams,mode=mode)
        self.set_priors()


    def set_logLevel(self,logLevel):
        self.log.setLevel(logLevel)


    def prepare_data(self,modelParams,mode='rates',key=None):
        '''
            preprocesses the data, provided as an nAnimal x max(nNeuron) array,
            containing the spike counts of each neuron

            nNeuron might differ between animals, so the data is padded with NaNs
        '''

        self.T = modelParams.T
        

        if mode=='rates':
            if key:
                self.rates = np.array(modelParams.rates[key])
            else:
                self.rates = np.array(modelParams.rates)
            self.spike_counts = np.floor(self.rates * self.T)
        else:
            if key:
                self.spike_counts = np.array(modelParams.spike_counts[key])
            else:
                self.spike_counts = np.array(modelParams.spike_counts)
            self.rates = self.spike_counts / self.T
        ## get the maximum number of spikes
        self.N_max = np.nanmax(self.spike_counts).astype('int64')
        self.nSamples = self.spike_counts.shape[1]      ## here: number of animals
        
        self.nAnimals = self.nSamples

        self.data = {
            # 'nAnimals' : self.spike_counts.shape[1],
            'nNeurons' : [np.isfinite(spikes).sum() for spikes in self.spike_counts.T],
        }

        ## calculate the number of occurences 'k_AP' of each spike number per animal,
        ## as well as the binomial coefficient of occurence 'k_AP', given n neurons
        ## can be sped up, using: 
        ##      if k=0, then the binomial coefficient is 1
        ##      if k=1, then the binomial coefficient is n
        self.data['k_AP'] = np.zeros(
            (self.nSamples,self.N_max+1)
        )
        self.binom = np.zeros_like(self.data['k_AP'])
        for a,N in enumerate(self.spike_counts.T):
        
            N = N[np.isfinite(N)].astype('int64')
            N_ct = Counter(N)
            self.data['k_AP'][a,list(N_ct.keys())] = list(N_ct.values())
        
            self.binom[a,:] = sp_binom(self.data['nNeurons'][a],self.data['k_AP'][a])


    def flatten_parameters(self, params):

        return params.reshape((params.shape[0],-1),order='F')
    

    def unflatten_parameters(self, params):

        return params.reshape((params.shape[0],-1,self.nParams),order='F')


    # def un_flatten_parameters(self,params, force_mode=None):
        
    #     N = params.shape[0]
    #     if len(params.shape)==1:
    #         return params[np.newaxis,...]
    #     elif len(params.shape)==2:
    #         if force_mode=='flatten':
    #             return params
    #         else:
    #             # print('unflatten:',params.shape)
    #             # assert params.shape[1]/self.nParams == self.data['nAnimals'], 'provided parameter array does not match the number of animals'
    #             return params.reshape((N,-1,self.nParams),order='C')
    #     elif len(params.shape)==3:
    #         if force_mode=='unflatten':
    #             return params
    #         else:
    #             # assert np.prod(params.shape[1:])/self.nParams == self.data['nAnimals'], 'provided parameter array does not match the number of animals'
    #             # print('flatten:',params.shape)
    #             return params.reshape((N,-1),order='C')
    #     else:
    #         return params


    def set_priors(self,hierarchical=[],two_pop=False):

        '''
            builds prior distributions for the model parameters
            each parameter can either be hierarchical or not
        '''

        self.two_pop = two_pop

        ## define distributions inline, instead of using scipy ones for the sake of performance
        halfnorm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(x)
        norm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(2*x - 1)
        
        self.priors_init = {
            'gamma': {
                # 'hierarchical':     False,
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 2., 'scale': 0.5},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.1},
                    'function':     halfnorm_ppf,
                },
            },
            'delta': {
                # 'hierarchical':     False,
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 4, 'scale': 1},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.3},
                    'function':     halfnorm_ppf,
                },
            },
            'nu_max': {
                # 'hierarchical':     False,
                'hierarchical':     {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 30, 'scale': 5},
                    'function':     norm_ppf,
                },
                'params':           {'loc': 'mean', 'scale': 'sigma'},
                'sigma': {
                    'params':       {'loc': 0, 'scale': 2},
                    'function':     halfnorm_ppf,
                },
            },
        }
        
        if self.two_pop:
            
            self.priors_init['weight_dark'] = {
                'mean': {
                    'params':       {},
                    'function':     lambda x : 0.5,
                },
            }

            self.priors_init['gamma_dark'] = {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 1., 'scale': 0.3},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.1},
                    'function':     halfnorm_ppf,
                },
            }

            self.priors_init['delta_dark'] = {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 4, 'scale': 1},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.3},
                    'function':     halfnorm_ppf,
                },
            }

        '''
            translating the above specified priors into a format that can be used by the model

            TODO:
                * this part of the function can be moved to a general function, that can be used for any model

        '''
        self.paraNames = []
        self.priors = {}
        # self.pTC = {}
        
        
        ct = 0
        for key in self.priors_init:

            if key in hierarchical:
                # self.priors_init[key]['hierarchical']:
                
                ## add the mean and sigma parameters for the hierarchical prior
                for var in self.priors_init[key]['hierarchical']['params'].values():
                    
                    paraName = f"{key}_{var}"
                    # print(paraName)
                    self.priors[paraName] = {}
                
                    self.paraNames.append(paraName)
                    
                    self.priors[paraName]['idx'] = ct
                    self.priors[paraName]['n'] = 1
                    self.priors[paraName]['meta'] = True

                    self.priors[paraName]['transform'] = \
                        lambda x,params=self.priors_init[key][var]['params'],fun=self.priors_init[key][var]['function']: fun(x,**params)
                    
                    ct += 1

                ## then, add the actual parameters for the hierarchical prior
                paraName = f"{key}"
                self.priors[paraName] = {}
                # print(paraName)
                for i in range(self.nSamples):
                    self.paraNames.append(f'{key}_{i}')
                self.priors[paraName]['idx'] = ct
                # get indexes of hierarchical parameters for quick access later on
                self.priors[paraName]['idx_mean'] = self.priors[f"{key}_{self.priors_init[key]['hierarchical']['params']['loc']}"]['idx']
                self.priors[paraName]['idx_sigma'] = self.priors[f"{key}_{self.priors_init[key]['hierarchical']['params']['scale']}"]['idx']
                
                self.priors[paraName]['n'] = self.nSamples
                self.priors[paraName]['meta'] = False

                self.priors[paraName]['transform'] = \
                    lambda x,params,fun=self.priors_init[key]['hierarchical']['function']: fun(x,**params)

                ct += self.nSamples
            
            else:
                var = "mean"
                paraName = f"{key}"

                self.paraNames.append(paraName)
                # print(paraName)
                self.priors[paraName] = {}

                self.priors[paraName]['idx'] = ct
                self.priors[paraName]['n'] = 1
                self.priors[paraName]['meta'] = False

                self.priors[paraName]['transform'] = \
                    lambda x,params=self.priors_init[key][var]['params'],fun=self.priors_init[key][var]['function']: fun(x,**params)
                ct += 1
        
        self.nParams = len(self.paraNames)
        self.wrap = np.zeros(self.nParams).astype('bool')
        # for key in self.priors:
        #     print(key)
        #     key_root, key_var = key.split('_')
        #     if 'hierarchical' in self.priors_init[key_root]:
        #         self.priors['wrap'][self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']] = self.priors_init[key_root][key_var]['wrap']
        

    
    def set_prior_transform(self,vectorized=True,mode='vectorized'):
        '''
            sets the prior transform function for the model, 
            which transforms the random variables from the unit hypercube to the actual priors

            only takes as input the mode, which can be either of
            - 'vectorized': vectorized prior transform function
            - 'scalar': scalar prior transform function
            - 'tensor': tensor prior transform function
        '''

    # def prior_transform(self,p_in,vectorized=True):
        def prior_transform(p_in):

            # def transform_p(self,p_in,vectorized=True):
        
            if len(p_in.shape)==1:
                p_in = p_in[np.newaxis,...]
            p_out = np.zeros_like(p_in)
            
            for key in self.priors:
                
                if self.priors[key]['n']==1:
                
                    p_out[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']] = self.priors[key]['transform'](p_in[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']])

                else:
                    params = {
                        'loc':          p_out[:,self.priors[key]['idx_mean'],np.newaxis],
                        'scale':        p_out[:,self.priors[key]['idx_sigma'],np.newaxis],
                        # 'loc':          0,
                        # 'scale':        1,
                    }
                    # print(key,self.priors[key])
                    # print(params)
                    # print(p_in[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']])
                    p_out[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']] = self.priors[key]['transform'](p_in[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']],params=params)
                    # print(p_out[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']])
            
            if vectorized:
                return p_out
            else:
                return p_out[0,:]
            
        return prior_transform
        

    def set_logl(self,vectorized=True,N_low_counts=1):
        '''
            sets the log likelihood function for the model,
            which calculates the log likelihood of the data given the model parameters

            it assumes that neuron activity is provided as spike counts per neuron
            and the model is a Poisson model
        '''


        def logl_count_based(p,N_AP_low,k_AP,binom,nNeurons,zero=False):

            '''
                calculates probability to observe k_AP neurons with spike count N_AP, 
                when drawing nNeurons from the distribution
                
                mostly used for low counts (for both observed, and non-observed counts)
                
                calculates
                    p_N_AP [nk]
                        The probability to observe N_AP action potentials in any neuron, given the underlying firing rate distribution and poisson firing. Calculated as the integral over all possible rates
                            int_0^{nu^max} p(nu | Sigma) * p(N_AP | nu,T) dnu
                    
                    p_k_AP [nk]
                        The probability to observe the empirical spike count distribution, given p_N_AP and the number of neurons. Calculated with the binomial distribution
                            binomial(n k) * p_N_AP^k * (1-p_N_AP)^(n-k)

                calculates as log, to avoid numerical issues (p->0, logp->-inf)
            '''

            # if zero:
                # p_N_AP = adaptive_integration(f,0,p['nu_max'],
                    # args=(*p.values(),0,self.T,True))
            # else:
            p_N_AP = adaptive_integration(f,0,p['nu_max'],
                args=(p,N_AP_low,self.T,zero,self.two_pop)
            )

            if p_N_AP is None:
                self.log.warning('error in integration for p_N_AP')
                return -100000.
            
            logp = np.log(binom) + k_AP * np.log(p_N_AP) + (nNeurons - k_AP) * np.log(1-p_N_AP)
            logp -= np.log(nNeurons * p['nu_max'])    ## normalization to density

            # print(f"{k_AP=}, {p_N_AP=} {logp_k=}, {logl_low=}")
            return logp.sum()
        

        def logl_continuous(p,rates):
            '''
                calculates the overall probability by evaluating empirical rates
                at theoretical probability density function to obtain probability of
                observing rates, given current model parameters
            '''
            p_continuous = p_nu(rates,p,two_pop=self.two_pop)
            logp = np.log(p_continuous) - np.log(p['nu_max']*self.T)        ## probability to observe firing range in small interval, defined by integer counts
            
            return logp.sum()


        def logl_tail(p,rates):
            
            '''
                finally, estimate probability, that no "more extreme" firing rates are observed when drawing nNeurons neurons from the distribution by calculating the integral from 0 to the highest observed firing rate
            '''            
            emp_nu_max = rates.max()
            nNeurons = len(rates)
            p_spike_count_in_empirical_range = adaptive_integration(p_nu,0,emp_nu_max,
                    args=(p,self.two_pop))
            
            if p_spike_count_in_empirical_range is None:
                self.log.warning('error in integration for tail')
                return -100000.

            logl_no_spikes_in_tail = nNeurons * np.log(p_spike_count_in_empirical_range)
            self.log.debug(f'up to {emp_nu_max=}: {p_spike_count_in_empirical_range=}, p_no_spikes_in_tail={np.exp(logl_no_spikes_in_tail)}')
            return np.maximum(logl_no_spikes_in_tail,-100000)
            


        def loglikelihood(params,withZeros=False,plot=False):
            '''
                log likelihood function for the model
                
                The likelihood is calculated as count based for low expected spike counts (nu_max*T < N_low_counts) 
                and as continuous for medium and high expected spike counts (nu_max*T >= N_low_counts)
                
                input (* = class-variable)
                    params [nChain, nAnimals, nParams]
                        The input parameters for the model, which are the gamma, delta and nu_max

                    * N_AP [nAnimals, nk]
                        The number of action potentials observed in each neuron
                    
                    * k_AP [nAnimals, nk]
                        The number of occurences of each spike count in the data per animal

                    * binom [nAnimals, nk]
                        The binomial coefficient to draw 'k' neurons with some property of probability p, given 'nNeuron' neurons

                output
                    logp [nChain]
                        The log likelihood of the data, given the model parameters
                            sum_{nAnimals} sum_{N_max} log( p_k_AP )

                dimensions:
                    nChain:                     
                        number of chains
                    nAnimals:
                        number of animals
                    nNeurons [nAnimals]:        
                        number of neurons (possibly different for each animal)
                    nParams:                    
                        number of parameters
                    nk [nChain, nAnimals]:   
                        maximum possible number of spikes given nu_max
            '''
            # define dimensions
            if len(params.shape)==1:
                params = params[np.newaxis,:]
            
            
            nChain = params.shape[0]
            
            logl = np.zeros((nChain,self.nAnimals))

            for a in range(self.nAnimals):
                
                if plot:
                    fig,ax = plt.subplots(2,nChain,figsize=(18,5))

                rates = self.rates[np.isfinite(self.rates[:,a]),a]
                spike_counts = self.spike_counts[np.isfinite(self.spike_counts[:,a]),a]


                # p = {}
                # for var in ['gamma','delta','nu_max']:
                #     p[var] = params[:,self.priors[var]['idx'] + ( a if self.priors[var]['n']>1 else 0 )][np.newaxis,...]

                for i in range(nChain):
                    
                    p = {}
                    for var in self.priors_init.keys():
                        p[var] = params[i,self.priors[var]['idx'] + ( a if self.priors[var]['n']>1 else 0 )]

                    if p['nu_max'] < rates.max():
                        self.log.warning('nu_max too low for data')
                        logl[i,a] = -10**6 * (rates.max() - p['nu_max'])
                        continue
                    if 'gamma_dark' in p and p['gamma_dark'] > p['gamma']:
                        self.log.warning('gamma_dark > gamma')
                        logl[i,a] = -10**6 * (p['gamma_dark'] - p['gamma'])
                        continue

                    self.log.debug(f"animal {a}, parameters: {p['gamma']=}, {p['delta']=}, {p['nu_max']=}")


                    nNeurons = self.data['nNeurons'][a]
                    
                    nk = int(min(N_low_counts,p['nu_max'] * self.T))

                    if nk > 0:
                        logl_low = logl_count_based(
                            p,np.arange(nk),
                            self.data['k_AP'][a,:nk],
                            self.binom[a,:nk],
                            nNeurons,
                            nk==1
                        )
                    else:
                        logl_low = 0.


                    '''
                        for other counts, use firing rate distribution (to save computational power)
                    '''
                    # N_AP = self.spike_counts[:,a]
                    # N_AP = N_AP[N_AP>=nk]
                    logl_intermediate = logl_continuous(p,rates=rates[spike_counts>nk])
                    

                    # N_AP = np.where(self.data['k_AP'][a,:]>0)[0]                    
                    # logl_tail = logl_tail(p,self.rates)
                    logl_tail = 0.
                    

                    '''
                        sum up the log likelihoods to obtain the overall estimate for this animal
                    '''
                    self.log.debug(('logls:',logl_low,logl_intermediate,logl_tail))

                    logl[i,a] = logl_low + logl_intermediate + logl_tail
                    if not np.isfinite(logl[i,a]):
                        logl[i,a] = -100000.

            self.log.debug(logl)

            if vectorized:
                return logl.sum(axis=1)
            else:
                return logl[0,:].sum()
        
        return loglikelihood



def adaptive_integration(f,x_lower,x_upper,args,eps_thr=-4):
    eps_pow = -8
    while True:
        if eps_pow==eps_thr:
            # self.log.warning('tolerance too high - breaking!')
            return None

        # try:
        res = quadpy.quad(f,x_lower,x_upper,args=args,epsabs=10**eps_pow, epsrel=10**eps_pow)[0]
        break
        # except:
            # self.log.debug(f'error in integration with tolerance 10^{eps_pow}')
            # eps_pow += 1
    return res


def run_sampling(mP,mode='ultranest',two_pop=False):

    BM = BayesModel(mP,mode='rates')
    BM.set_logLevel(logging.ERROR)
    BM.prepare_data(mP,mode='rates',key='WT')
    BM.set_priors(two_pop=two_pop)
    
    my_prior_transform = BM.set_prior_transform()
    my_likelihood = BM.set_logl(N_low_counts=1)


    if mode=='dynesty':
        my_prior_transform = BM.set_prior_transform(vectorized=False)
        my_likelihood = BM.set_logl(vectorized=False)
        # my_prior_transform = lambda p_in : hbm.transform_p(p_in,vectorized=False)
        # my_likelihood = hbm.set_logl_func(vectorized=False)
        print('running nested sampling')
        # print(np.where(hbm.pTC['wrap'])[0])
        with dypool.Pool(8,my_likelihood,my_prior_transform) as pool:
            sampler = DynamicNestedSampler(pool.loglike,pool.prior_transform,BM.nParams,
                    pool=pool,
                    # periodic=np.where(BM.wrap)[0],
                    sample='slice'
                )
            sampler.run_nested()

        sampling_result = sampler.results
        print(sampling_result)
        return BM, sampling_result, sampler
    else:
        # print(hbm.paraNames)
        sampler = ultranest.ReactiveNestedSampler(
            BM.paraNames, 
            my_likelihood, my_prior_transform,
            wrapped_params=BM.wrap,
            vectorized=True,num_bootstraps=20,
            ndraw_min=512
        )

        logger = logging.getLogger("ultranest")
        logger.setLevel(logging.ERROR)

        nsteps = 2*BM.nParams
        sampler.stepsampler = ultranest.stepsampler.SliceSampler(
            nsteps=nsteps,
            generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        )


        # num_samples = BM.nParams*100
        # num_samples = np.maximum(400,BM.nParams*100)
        num_samples = 400

        sampling_result = sampler.run(
            min_num_live_points=num_samples,
            max_iters=20000,cluster_num_live_points=20,max_num_improvement_loops=3,
            show_status=True,viz_callback='auto')
        
        # plt.figure()
        # cornerplot(sampling_result)
        # plt.show(block=False)

        # plt.figure()
        # sampler.plot_trace()
        # plt.show(block=False)

        return BM, sampling_result, sampler

