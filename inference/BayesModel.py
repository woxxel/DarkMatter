import os, logging
import numpy as np

from scipy.stats import truncnorm
from scipy.special import binom as sp_binom, erfinv
from scipy.special import gammaln#factorial as sp_factorial

from collections import Counter
import quadpy

from matplotlib import pyplot as plt

logging.basicConfig(level=logging.INFO)



def p_nu(NU,gamma,delta,nu_max,log=False):
    
    if log:
        return - np.log( nu_max / gamma * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) - delta**2 / 2 + \
        ( gamma**2 - 1 ) * np.log( NU / nu_max ) + \
        np.log( np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max ) ) ) )
    else:
        return gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) * \
            np.exp( - delta**2/2.) * ( NU / nu_max )**(gamma**2 - 1) * \
            np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max) ) )

def poisson_spikes(nu,N_AP,T,log=False):
    ## using the gamma-function to obtain log(N!) = gammaln(N+1)
    if log:
        return N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T
    else:
        return np.exp(N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T)


def f(nu,gamma,delta,nu_max,N_AP,T,loop=False,log=False):
    if loop:
        return p_nu(nu,gamma,delta,nu_max,log=log) * poisson_spikes(nu,N_AP,T,log=log)
    else:
        # print('poisson:',poisson_spikes(nu,N_AP[:,np.newaxis],T,log=log)[:200,:])
        return p_nu(nu,gamma,delta,nu_max,log=log) * poisson_spikes(nu,N_AP[:,np.newaxis],T,log=log)



class BayesModel:


    def __init__(self,modelParams,logLevel=logging.ERROR):
        '''
            initialize the model by setting the priors and preparing the data
        '''

        self.log = logging.getLogger("nestLogger")
        self.log.setLevel(logLevel)


        os.environ['MKL_NUM_THREADS'] = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        os.environ['OMP_NUM_THREADS'] = '1'

        self.param_names = ['gamma','delta','nu_max']
        
        self.prior = {
            'gamma': {'mu': 1.5, 'std':0.3, 'lower':0., 'upper': np.inf},
            'delta': {'mu': 4., 'std':1., 'lower': 0., 'upper': np.inf},
            'nu_max': {'mu': 35., 'std':2., 'lower': 5., 'upper': np.inf},
        }
        # self.nParams = len(self.prior)

        self.T = modelParams.T
        self.N = modelParams.N
        self.spike_counts = np.array(modelParams.spike_counts)
        
        self.prepare_data()


    def prepare_data(self):
        '''
            preprocesses the data, provided as an nAnimal x max(nNeuron) array,
            containing the spike counts of each neuron

            nNeuron might differ between animals, so the data is padded with NaNs
        '''

        ## get the maximum number of spikes
        self.N_max = np.nanmax(self.spike_counts).astype('int64')

        self.spike_counts = self.spike_counts
        
        self.nSamples = self.spike_counts.shape[1]      ## here: number of animals

        self.data = {
            'nAnimals' : self.spike_counts.shape[1],
            'nNeurons' : [np.isfinite(spikes).sum() for spikes in self.spike_counts.T],
        }

        ## calculate the number of occurences 'k_AP' of each spike number per animal,
        ## as well as the binomial coefficient of occurence 'k_AP', given n neurons
        ## can be sped up, using: 
        ##      if k=0, then the binomial coefficient is 1
        ##      if k=1, then the binomial coefficient is n
        self.data['k_AP'] = np.zeros(
            (self.data['nAnimals'],self.N_max+1)
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


    def set_priors(self):

        '''
            builds prior distributions for the model parameters
            each parameter can either be hierarchical or not
        '''

        ## define distributions inline, instead of using scipy ones for the sake of performance
        halfnorm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(x)
        norm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(2*x - 1)
        
        self.priors_init = {
            'gamma': {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 1.5, 'scale': 0.3},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.1},
                    'function':     halfnorm_ppf,
                },
            },
            'delta': {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 3, 'scale': 1},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.3},
                    'function':     halfnorm_ppf,
                },
            },
            'nu_max': {
                'hierarchical':     False,
                # 'hierarchical':     norm_ppf,
                'mean': {
                    'params':       {'loc': 30, 'scale': 5},
                    'function':     norm_ppf,
                },
                # 'params':           {'loc': 'mean', 'scale': 'sigma'},
                # 'sigma': {
                #     'params':       {'loc': 0, 'scale': 1},
                #     'function':     halfnorm_ppf,
                # },
            },
        }
        

        '''
            translating the above specified priors into a format that can be used by the model
        '''
        self.paraNames = []
        self.priors = {}
        # self.pTC = {}
        
        
        ct = 0
        for key in self.priors_init:

            if self.priors_init[key]['hierarchical']:
                
                ## add the mean and sigma parameters for the hierarchical prior
                for var in self.priors_init[key]['hierarchical']['params'].values():
                    
                    paraName = f"{key}_{var}"
                    print(paraName)
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
                print(paraName)
                for i in range(self.nSamples):
                    self.paraNames.append(f'{key}_{i}')
                self.priors[paraName]['idx'] = ct
                # get indexes of hierarchical parameters for quick access later on
                self.priors[paraName]['idx_mean'] = self.priors[f"{key}_{self.priors_init[key]['hierarchical']['params']['loc']}"]['idx']
                self.priors[paraName]['idx_sigma'] = self.priors[f"{key}_{self.priors_init[key]['hierarchical']['params']['scale']}"]['idx']
                
                self.priors[paraName]['n'] = self.nSamples
                self.priors[paraName]['meta'] = False

                self.priors[paraName]['transform'] = \
                    lambda x,params,fun=self.priors_init[key][var]['function']: fun(x,**params)

                ct += self.nSamples
            
            else:
                var = "mean"
                paraName = f"{key}"

                self.paraNames.append(paraName)
                print(paraName)
                self.priors[paraName] = {}

                self.priors[paraName]['idx'] = ct
                self.priors[paraName]['n'] = 1
                self.priors[paraName]['meta'] = False

                self.priors[paraName]['transform'] = \
                    lambda x,params=self.priors_init[key][var]['params'],fun=self.priors_init[key][var]['function']: fun(x,**params)
        
        self.nParams = len(self.paraNames)
        self.wrap = np.zeros(self.nParams).astype('bool')
        # for key in self.priors:
        #     print(key)
        #     key_root, key_var = key.split('_')
        #     if 'hierarchical' in self.priors_init[key_root]:
        #         self.priors['wrap'][self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']] = self.priors_init[key_root][key_var]['wrap']
        

    
    def set_prior_transform(self,mode='vectorized'):
        '''
            sets the prior transform function for the model, 
            which transforms the random variables from the unit hypercube to the actual priors

            only takes as input the mode, which can be either of
            - 'vectorized': vectorized prior transform function
            - 'scalar': scalar prior transform function
            - 'tensor': tensor prior transform function
        '''

    # def prior_transform(self,p_in,vectorized=True):
        def prior_transform(p_in,vectorized=True):

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
                    p_out[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']] = self.priors[key]['transform'](p_in[:,self.priors[key]['idx']:self.priors[key]['idx']+self.priors[key]['n']],params=params)
            
            if vectorized:
                return p_out
            else:
                return p_out[0,:]
            
        return prior_transform

        # def prior_transform(cube):

        #     '''
        #         receives N x (nParams x nAnimals) parameter-array in unit hypercube space
        #         to transform into the actual prior space

        #         g,g,g,g,g,d,d,d,d,d,n,n,n,n,n
        #     '''

        #     N = cube.shape[0]
        #     params_in = cube.reshape(-1,self.nParams)
            
        #     print(params_in)
        #     # transforms random variables from uniform [0,1] distribution into actual priors
        #     params_out = np.zeros_like(params_in)
        #     for i,var in enumerate(self.prior.keys()):
        #         mu = self.prior[var]['mu']
        #         sigma = self.prior[var]['std']
        #         lower = self.prior[var]['lower']
        #         upper = self.prior[var]['upper']
        #         params_out[...,i] = truncnorm.ppf(params_in[...,i],(lower-mu)/sigma,(upper-mu)/sigma,mu,sigma)
            
        #     return params_out.reshape(N,)
        # return prior_transform
        

    def set_logl(self,mode='vectorized',N_low_counts=100):
        '''
            sets the log likelihood function for the model,
            which calculates the log likelihood of the data given the model parameters

            it assumes that neuron activity is provided as spike counts per neuron
            and the model is a Poisson model
        '''

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
            if len(params.shape)==2:
                params = params[np.newaxis,:]
            
            nChain = params.shape[0]
            nAnimals = self.data['nAnimals']

            logl = np.zeros((nChain,nAnimals))

            for i,params_ in enumerate(params):
                
                if plot:
                    fig,ax = plt.subplots(2,nAnimals,figsize=(18,5))
                for a, (gamma,delta,nu_max) in enumerate(params_):
                    # print('parameters: ',gamma,delta,nu_max)
                    if (nu_max*self.T) < self.spike_counts[a,:].max():
                        continue

                    nNeurons = self.data['nNeurons'][a]
                    ## allocate some resources
                    

                    '''
                        for low counts, use count-based inference (for both observed, and non-observed counts)
                        
                        calculates
                            p_N_AP [nk]
                                The probability to observe N_AP action potentials in any neuron, given the underlying firing rate distribution and poisson firing. Calculated as the integral over all possible rates
                                    int_0^{nu^max} p(nu | Sigma) * p(N_AP | nu,T) dnu
                            
                            p_k_AP [nk]
                                The probability to observe the empirical spike count distribution, given p_N_AP and the number of neurons. Calculated with the binomial distribution
                                    binomial(n k) * p_N_AP^k * (1-p_N_AP)^(n-k)
                    '''
                    nk = int(min(N_low_counts,nu_max * self.T))
                    if N_low_counts > 0:
                        N_AP_low = np.arange(nk)
                        binom = self.binom[a,:nk]
                        k_AP = self.data['k_AP'][a,:nk]
                        

                        p_N_AP = adaptive_integration(f,0,nu_max,
                                args=(gamma,delta,nu_max,N_AP_low,self.T,False))

                        if p_N_AP is None:
                            self.log.warning('error in integration for p_N_AP')
                            logl[i,a] = -10000
                            continue
                        
                        ## calculate probability to observe k_AP neurons with spike count N_AP, when drawing nNeurons from the distribution
                        p_k = binom * p_N_AP**k_AP * (1-p_N_AP)**(nNeurons - k_AP)
                        p_k /= (nNeurons * nu_max)
                        logl_low = np.log(p_k).sum()
                    else:
                        logl_low = 0

                    '''
                        for other counts, use firing rate distribution (to save computational power)
                    '''
                    N_AP = np.where(self.data['k_AP'][a,:]>0)[0]
                    p_continuous = p_nu(N_AP/self.T,gamma,delta,nu_max)
                    p_steps = p_continuous / (nu_max*self.T)        ## probability to observe firing range in small interval, defined by integer counts
                    logl_intermediate = np.log(p_steps[N_AP>=nk]).sum()
                    

                    '''
                        finally, estimate probability, that no "more extreme" firing rates are observed when drawing nNeurons neurons from the distribution by calculating the integral from 0 to the highest observed firing rate
                    '''
                    emp_nu_max = self.spike_counts[a,:].max()/self.T
                    p_spike_count_in_empirical_range = adaptive_integration(p_nu,0,emp_nu_max,
                            args=(gamma,delta,nu_max,False))
                    if p_spike_count_in_empirical_range is None:
                        self.log.warning('error in integration for tail')
                        logl[i,a] = -10000
                        continue

                    p_no_spikes_in_tail = p_spike_count_in_empirical_range**nNeurons
                    self.log.debug(f'parameters: {gamma=}, {delta=}, {nu_max=}')
                    self.log.debug(f'up to {emp_nu_max=}: {p_spike_count_in_empirical_range=}, {p_no_spikes_in_tail=}')
                    logl_tail = np.log(p_no_spikes_in_tail)
                    

                    '''
                        sum up the log likelihoods to obtain the overall estimate for this animal
                    '''
                    self.log.debug(('logls:',logl_low,logl_intermediate,logl_tail))

                    logl[i,a] = logl_low + logl_intermediate + logl_tail


                    if plot:
                        line, = ax[0][a].plot(N_AP/self.T,p_steps)
                        col = line.get_color()
                        
                        ax[1][a].plot(N_AP/self.T,p_steps,'x',markersize=4,color=col)
                        if N_low_counts > 0:
                            ax[1][a].plot(N_AP_low[k_AP==0]/self.T,p_k[k_AP==0],'o',markersize=2,color='k')
                            for k in range(1,5):
                                ax[1][a].plot(N_AP_low[k_AP==k]/self.T,p_k[k_AP==k],'o',markersize=2,color=[1-k/5,0,0])

                        plt.setp(ax[0][a],ylim=(0,0.0002))
                        plt.setp(ax[1][a],ylim=(0,0.0002))

                if plot:
                    plt.show(block=False)
            self.log.debug(logl)
            # logl = np.zeros(nChain)
            return logl.sum(axis=1)
        return loglikelihood



def adaptive_integration(f,x_lower,x_upper,args,eps_thr=-4):
    eps_pow = -8
    while True:
        if eps_pow==eps_thr:
            # self.log.warning('tolerance too high - breaking!')
            return None

        try:
            res = quadpy.quad(f,x_lower,x_upper,args=args,epsabs=10**eps_pow, epsrel=10**eps_pow)[0]
            break
        except:
            # self.log.debug(f'error in integration with tolerance 10^{eps_pow}')
            eps_pow += 1
    return res