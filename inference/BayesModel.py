import os, logging, time
import numpy as np

from scipy.stats import truncnorm
from scipy.special import binom as sp_binom, erfinv
from scipy.special import gammaln#factorial as sp_factorial
from scipy.integrate import quad

from collections import Counter
import quadpy

from matplotlib import pyplot as plt

import ultranest
import ultranest.stepsampler
from ultranest.plot import cornerplot
from dynesty import DynamicNestedSampler, pool as dypool, utils, plotting as dyplot

from .utils.utils import adaptive_integration, p_nu, f

from DM_theory.functions import get_nu_bar, get_q

import warnings
warnings.filterwarnings("ignore")


logging.basicConfig(level=logging.INFO)



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
        # self.binom = np.zeros_like(self.data['k_AP'])
        self.binom = {}
        for a,N in enumerate(self.spike_counts.T):
        
            N = N[np.isfinite(N)].astype('int64')
            N_ct = Counter(N)
            self.data['k_AP'][a,list(N_ct.keys())] = list(N_ct.values())
            

            self.binom[a] = {}
            self.binom[a]['N_AP'] = np.where(self.data['k_AP'][a]>0)[0]
            self.binom[a]['k_AP'] = self.data['k_AP'][a][self.binom[a]['N_AP']].astype('int')
            self.binom[a]['binomial_coefficient'] = sp_binom(self.data['nNeurons'][a],np.arange(self.data['nNeurons'][a]+1))

            # self.binom[a] = binom[self.data['k_AP'][a]]

        ## predefine for calculations in logl function
        # self.k_idxes = np.where(self.data['k_AP']>0)
        # self.k_raw = self.data['k_AP'][self.k_idxes]



    # def flatten_parameters(self, params):
    #     return params.reshape((params.shape[0],-1),order='F')
    

    # def unflatten_parameters(self, params):
    #     return params.reshape((params.shape[0],-1,self.nParams),order='F')


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
                    'params':       {'loc': 6, 'scale': 2},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.5},
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
            
            self.priors_init['p'] = {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },
                
                'mean': {
                    'params':       {},
                    'function':     lambda x : x,
                },
                'sigma': {
                    'params':       {'loc': 0., 'scale': 0.1},
                    'function':     norm_ppf,
                },
            }

            self.priors_init['gamma2'] = {
            # self.priors_init['delta_dark'] = {
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
            }

            self.priors_init['delta2'] = {
                'hierarchical': {
                    'function':     norm_ppf,
                    'params':       {'loc': 'mean', 'scale': 'sigma'},
                },

                'mean': {
                    'params':       {'loc': 6, 'scale': 2},
                    'function':     norm_ppf,
                },
                'sigma': {
                    'params':       {'loc': 0, 'scale': 0.5},
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
        

    def set_logl(self,vectorized=False,withZeros=False,N_low_counts=1):
        '''
            sets the log likelihood function for the model,
            which calculates the log likelihood of the data given the model parameters

            it assumes that neuron activity is provided as spike counts per neuron
            and the model is a Poisson model
        '''


        # def logl_count_based(p,N_AP_low,k_AP,binom,nNeurons,zero=False):

        #     '''
        #         calculates probability to observe k_AP neurons with spike count N_AP, 
        #         when drawing nNeurons from the distribution
                
        #         mostly used for low counts (for both observed, and non-observed counts)
                
        #         calculates
        #             p_N_AP [nk]
        #                 The probability to observe N_AP action potentials in any neuron, given the underlying firing rate distribution and poisson firing. Calculated as the integral over all possible rates
        #                     int_0^{nu^max} p(nu | Sigma) * p(N_AP | nu,T) dnu
                    
        #             p_k_AP [nk]
        #                 The probability to observe the empirical spike count distribution, given p_N_AP and the number of neurons. Calculated with the binomial distribution
        #                     binomial(n k) * p_N_AP^k * (1-p_N_AP)^(n-k)

        #         calculates as log, to avoid numerical issues (p->0, logp->-inf)
        #     '''

        #     zero = False
        #     # if zero:
        #         # p_N_AP = adaptive_integration(f,0,p['nu_max'],
        #             # args=(*p.values(),0,self.T,True))
        #     # else:
        #     p_N_AP = adaptive_integration(f,0,p['nu_max'],
        #         args=(p,N_AP_low,self.T,zero,self.two_pop)
        #     )

        #     # p_normalized = adaptive_integration(p_nu,0,p['nu_max'],
        #     #     args=(p,False)
        #     # )

        #     # print(N_AP_low)
        #     # p_continuous = p_nu(N_AP_low/self.T,p,two_pop=self.two_pop)

        #     # print(f"{p_N_AP=}, {p_continuous=}")#, {p_normalized=}")

        #     if p_N_AP is None:
        #         self.log.warning('error in integration for p_N_AP')
        #         return -100000.
            
        #     # p_N_AP /= (nNeurons * p['nu_max'])
        #     # print(np.log(nNeurons * p['nu_max']))
        #     logp = np.log(binom) + k_AP * np.log(p_N_AP) + (nNeurons - k_AP) * np.log(1-p_N_AP)
        #     # logp -= np.log(nNeurons * p['nu_max'])    ## normalization to density

        #     self.log.debug(f"count_based: {k_AP=}, {p_N_AP=}, {logp=}")
        #     # print(f"{k_AP=}, {p_N_AP=} {logp_k=}, {logl_low=}")
        #     return logp.sum()
        

        # def logl_continuous(p,rates):
        #     '''
        #         calculates the overall probability by evaluating empirical rates
        #         at theoretical probability density function to obtain probability of
        #         observing rates, given current model parameters
        #     '''
        #     p_continuous = p_nu(rates,p,two_pop=self.two_pop)
        #     logp = np.log(p_continuous) - np.log(self.T)        ## probability to observe firing range in small interval, defined by integer counts
        #     # self.log.debug(f"continuous (rates): {rates=}")
        #     # self.log.debug(f"continuous (logp): {logp=}")
        #     return logp.sum()


        # def logl_tail(p,rates):
            
        #     '''
        #         finally, estimate probability, that no "more extreme" firing rates are observed when drawing nNeurons neurons from the distribution by calculating the integral from 0 to the highest observed firing rate
        #     '''            
        #     emp_nu_max = rates.max()
        #     nNeurons = len(rates)
        #     p_spike_count_in_empirical_range = adaptive_integration(p_nu,0,emp_nu_max,
        #             args=(p,self.two_pop))
            
        #     if p_spike_count_in_empirical_range is None:
        #         self.log.warning('error in integration for tail')
        #         return -100000.

        #     logl_no_spikes_in_tail = nNeurons * np.log(p_spike_count_in_empirical_range)
        #     self.log.debug(f'up to {emp_nu_max=}: {p_spike_count_in_empirical_range=}, p_no_spikes_in_tail={np.exp(logl_no_spikes_in_tail)}')
        #     return np.maximum(logl_no_spikes_in_tail,-100000)
            


        def loglikelihood(params,plot=False):
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
            
            

            ## should be done somewhere else
            max_spike_count = np.nanmax(self.spike_counts,axis=0).astype('int')
            # print(max_spike_count)



            nChain = params.shape[0]
            
            logl = np.zeros((nChain,self.nAnimals))

            # k_idxes = np.where(self.data['k_AP']>0)
            # k_raw = self.data['k_AP'][k_idxes]

            for a in range(self.nAnimals):
                
                if plot:
                    fig,ax = plt.subplots(2,nChain,figsize=(18,5))

                rates = self.rates[np.isfinite(self.rates[:,a]),a]

                nNeurons = self.data['nNeurons'][a]
                max_spike_count = int(np.nanmax(self.spike_counts[:,a]))

                for i in range(nChain):
                    
                    p = {}
                    for var in self.priors_init.keys():
                        p[var] = params[i,self.priors[var]['idx'] + ( a if self.priors[var]['n']>1 else 0 )]
                    
                    # check_and_adjust_parameters(p)
                    # def check_and_adjust_parameters(p):

                    #     logl_penalty = 0.
                    if self.two_pop:
                        p['nu_max2'] = p['nu_max']

                        if p['gamma'] < p['gamma2']:
                            logl[i,a] += -10**6 * (p['gamma2'] - p['gamma'])
                        # p['weight_dark'] = 0.5

                        # p['delta_dark'] = np.sqrt(- (1+p['gamma_dark']**2) * np.log(p['gamma']**2/p['gamma_dark']**2) - (1 + p['gamma_dark']**2) * np.log((1+p['gamma_dark']**2)/(1+p['gamma']**2)) + (1+p['gamma_dark']**2)/(1+p['gamma']**2) * p['delta']**2)
                        # p['delta_dark'] = p['gamma_dark']/p['gamma'] * p['delta']
                        # p['gamma_dark'] = p['gamma']

                        ### enter constraints on parameters, especially gamma & delta!

                        if np.isnan(p['delta2']):
                            self.log.warning('delta2 is NaN')
                            logl[i,a] += -10**6
                            continue
                        
                        # return p, logl_penalty
                    

                    N_max = int(p['nu_max'] * self.T)

                    offset = 0.5
                    N_AP_array = np.arange(N_max).astype('int')
                    
                    N_AP_empirical = self.binom[a]['N_AP'][self.binom[a]['N_AP']<N_max]
                    k_AP_empirical = self.binom[a]['k_AP'][self.binom[a]['N_AP']<N_max]
                    
                    if withZeros:
                        N_AP = N_AP_array

                        k_AP = np.zeros(N_max,dtype='int')
                        k_AP[N_AP_empirical] = k_AP_empirical
                        
                    else:
                        N_AP = N_AP_empirical
                        k_AP = k_AP_empirical
                    
                    ## get all spike counts higher than allowed by model
                    idx_too_high = np.where(self.binom[a]['N_AP'] > N_max)[0]
                    for idx in idx_too_high:
                        # print(idx)
                        N_AP_too_high = self.binom[a]['N_AP'][idx]
                        k_AP_too_high = self.binom[a]['k_AP'][idx]
                        # print(f'{N_AP_too_high=}, {k_AP_too_high=}') 
                        logl[i,a] += -10**2 * k_AP_too_high*(N_AP_too_high - N_max)**2

                    # logl[i,a] += -100*((rates.mean() - get_nu_bar(p['gamma'],p['delta'],p['nu_max']))/rates.mean())**2
                    

                    # if 'gamma_dark' in p and p['gamma_dark'] > p['gamma']:
                    #     self.log.warning('gamma_dark > gamma')
                    #     logl[i,a] = -10**6 * (p['gamma_dark'] - p['gamma'])
                    #     continue

                    # self.log.debug(f"animal {a}, parameters: {p['gamma']=}, {p['delta']=}, {p['nu_max']=}")
                    

                    ## calculate maximum number of spikes:
                    binom = self.binom[a]['binomial_coefficient'][k_AP]

                    p_N_AP = p_nu((N_AP_array+offset)/self.T,p,two_pop=self.two_pop) / self.T

                    logp = np.log(binom) + k_AP * np.log(p_N_AP[N_AP]) + (nNeurons - k_AP) * np.log(1-p_N_AP[N_AP])

                    # print(f'{k_AP_empirical=}, {logp[N_AP_empirical]=}')

                    logl_measures = logp.sum()
                    logl[i,a] += logl_measures
                    

                    bias_to_mean = True
                    if bias_to_mean:
                        # adjust to two-population model!


                        if self.two_pop:
                            nu_mean = p['p'] * get_nu_bar(p['gamma'],p['delta'],p['nu_max'])
                            nu_mean += (1-p['p']) * get_nu_bar(p['gamma2'],p['delta2'],p['nu_max2'])
                            nu_SD = p['p'] * get_q(p['gamma'],p['delta'],p['nu_max'])
                            nu_SD += (1-p['p']) * get_q(p['gamma2'],p['delta2'],p['nu_max2'])
                        else:
                            nu_mean = get_nu_bar(p['gamma'],p['delta'],p['nu_max'])
                            nu_SD = get_q(p['gamma'],p['delta'],p['nu_max'])
                         
                        nu_sigma = np.sqrt(nu_SD - nu_mean**2) / np.sqrt(nNeurons)
                        # print(f'{nu_mean=}, {nu_SD=}, {nu_sigma=}')

                        logl_bias_to_mean = -(rates.mean() - nu_mean)**2 / (2*nu_sigma**2)
                        logl[i,a] += 10*logl_bias_to_mean
                    else:
                        logl_bias_to_mean = 0



                    extrema = True
                    if extrema and (max_spike_count < N_max):
                        ### now, get extreme value distribution
                        p_N_AP_cum = np.pad(
                            np.cumsum(p_N_AP)**nNeurons,
                            (1,0),mode='constant',constant_values=0
                        )
                        p_extreme = np.diff(p_N_AP_cum)
                        p_extreme *= (1 - 1./(p['nu_max']-N_AP_array/self.T+1)**2)
                        
                        logp_extreme_empirical = np.log(p_extreme[max_spike_count])
                        logl[i,a] += 10*logp_extreme_empirical
                    else:
                        logp_extreme_empirical = -10**2
                        logl[i,a] += logp_extreme_empirical

                    # print(f'{logl_measures=}, {logl_bias_to_mean=}, {logp_extreme_empirical=}')
                    # print(f'{logl[i,a]=}')

                    '''
                        sum up the log likelihoods to obtain the overall estimate for this animal
                    '''
                    # self.log.debug((f'logls: {logl_low=}, {logl_intermediate=}, {logl_tail=}'))

                    # logl[i,a] = logl_low + logl_intermediate + logl_tail
                    if not np.isfinite(logl[i,a]):
                        logl[i,a] = -100000.

            self.log.debug(logl)

            
            # if vectorized:
                # logl_sum = logl.sum(axis=1)
            return logl.sum(axis=1)
            # else:
                # return logl[0,:].sum()
        
        return loglikelihood



def run_sampling(mP,mode='ultranest',two_pop=False,withZeros=True):

    BM = BayesModel(mP,mode='rates')
    BM.set_logLevel(logging.ERROR)
    # BM.prepare_data(mP,mode='rates',key='WT')
    BM.prepare_data(mP,mode='rates')
    BM.set_priors(hierarchical=[],two_pop=two_pop)
    # BM.set_priors(hierarchical=['gamma','delta'],two_pop=two_pop)
    
    my_prior_transform = BM.set_prior_transform()
    my_likelihood = BM.set_logl(withZeros=withZeros)


    if mode=='dynesty':
        my_prior_transform = BM.set_prior_transform(vectorized=False)
        my_likelihood = BM.set_logl(vectorized=False,withZeros=withZeros)
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

        # nsteps = 2*BM.nParams
        # sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        #     nsteps=nsteps,
        #     generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        # )


        # num_samples = BM.nParams*100
        # num_samples = np.maximum(400,BM.nParams*100)
        num_samples = 1000

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



