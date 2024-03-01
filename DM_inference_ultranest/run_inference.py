import os, sys
root_dir = os.path.dirname(os.path.abspath(''))
if not root_dir in sys.path: sys.path.append(root_dir)

import numpy as np
import quadpy
import ultranest
from ultranest.plot import cornerplot

import matplotlib.pyplot as plt
from collections import Counter

import scipy as sp
from scipy.special import binom as sp_binom
from scipy.special import factorial as sp_factorial

import time

from inference import *
from DM_theory import *
from empirical.readData import *
from empirical.model import *

def p_nu(nu,gamma,delta,nu_max):

    #print(f"gamma={gamma}, delta={delta}, nu_max={nu_max}")
    return gamma / ( nu_max * np.sqrt( -np.pi * np.log( nu / nu_max ) ) ) * \
        np.exp( - delta**2/2.) * ( nu / nu_max )**(gamma**2 - 1) * \
        np.cosh( gamma * delta * np.sqrt( -2 * np.log( nu / nu_max) ) )

def poisson_spikes(nu,N,T_total):
    #print("poisson:",nu,T_total)
    return np.exp(N*np.log(nu*T_total) - np.log(sp_factorial(N)) - nu*T_total)
    #return (nu*T_total)**N / sp_factorial(N) * np.exp(-nu*T_total)

class Inference:

    dataLoaded = False
    paras = {}

    def __init__(self):

        os.environ['MKL_NUM_THREADS'] = '4'
        os.environ['OPENBLAS_NUM_THREADS'] = '4'
        os.environ['OMP_NUM_THREADS'] = '4'

        self.param_names = ['gamma','delta','nu_max']

        self.prior = {
            'gamma': {'mu': 1.5, 'std':0.5, 'lower':0., 'upper': np.inf},
            'delta': {'mu': 4, 'std':0.5, 'lower': 0., 'upper': np.inf},
            'nu_max': {'mu': 30., 'std':20, 'lower': 1., 'upper': np.inf},
        }
        self.ct=0


    def load_data(self, dataType='empirical', filePath='../data/BuscheLab/2P_data.xlsx',include_silent=False):

        self.mP = ModelParams(dataType, filePath=filePath, population_keys=['*mouse_type','animal'])
        self.data_df = self.mP.regularize_rates()

        self.data = self.data_df.to_numpy()

        N_zeros = (self.data==0).sum()
        #print(f'zeros in data: {N_zeros}')
#        self.data[self.data==0] = np.random.rand(N_zeros)*1./600
        #if include_silent:
        #    T = 600.
        #    self.data[self.data<=1./T] = -np.log(1-np.random.rand((self.data<=1./T).sum()))/T

        self.data_mask = ~np.isnan(self.data) & (self.data>0)

        self.dataLoaded = True


        ## prepare data such that we have all necessary values
        data_observed = self.prepare_data(0)


    def prepare_data(self,animal=0,T_total=10.,withZeros=True):

        N_AP = self.data*T_total

        N_AP_max = int(np.nanmax(N_AP))+1
        data_observed = {'n_animals': [],
                         #'N_AP': np.zeros((N_AP_max,N_AP.shape[1]),dtype="int64"),
                         'k_animals': np.zeros((N_AP.shape[1],N_AP_max),dtype="int64")
                        }

        for a,N in enumerate(N_AP.T):

            N = N[~np.isnan(N)].astype('int64')
            N_ct = Counter(N)
            data_observed['k_animals'][a,list(N_ct.keys())] = list(N_ct.values())

            data_observed['n_animals'].append(len(N))

        #rint(f'data: {self.data}')
        #rint(f'data observed: {data_observed}')
        #self.data_obs = data_observed

        self.n = data_observed['n_animals'][animal]

        if withZeros:
            self.k_AP = data_observed['k_animals'][animal]#[self.N_AP]#,np.newaxis]
            self.N_AP = np.arange(0,N_AP_max)[:,np.newaxis]
        else:
            self.N_AP = np.where(data_observed['k_animals'][animal])[0]#[:,np.newaxis]
            self.k_AP = data_observed['k_animals'][animal,self.N_AP]#,np.newaxis]
            self.N_AP = self.N_AP[:,np.newaxis]

        #print(self.N_AP.shape,self.k_AP.shape)
        self.T_total = T_total

        return data_observed

    def setPriorTransform(self):
        def prior_transform(cube):
            # transforms random variables from uniform [0,1] distribution into actual priors


            params = np.zeros_like(cube)
            for i,var in enumerate(self.prior.keys()):
                mu = self.prior[var]['mu']
                sigma = self.prior[var]['std']
                lower = self.prior[var]['lower']
                upper = self.prior[var]['upper']
                params[...,i] = sp.stats.truncnorm.ppf(cube[...,i],(lower-mu)/sigma,(upper-mu)/sigma,mu,sigma)
            #delta = sp.stats.truncnorm.ppf(cube[1],0,np.inf,self.prior['delta']['mu'],self.prior['delta']['std'])
            #nu_max = sp.stats.truncnorm.ppf(cube[2],0,np.inf,self.prior['nu_max']['mu'],self.prior['nu_max']['std'])

            #return [gamma,delta,nu_max]
            return params
        return prior_transform

    def setLogLikelihood(self,loop=True):
        def loglikelihood(params):
            ## define likelihood based on parameters, with data provided from context

            def f(nu,gamma,delta,nu_max,N_AP,T_total):
                #print(nu_max)
                if loop:
                    return p_nu(nu,gamma,delta,nu_max) * poisson_spikes(nu,N_AP,T_total)
                else:
                    return p_nu(nu,gamma,delta,nu_max) * poisson_spikes(nu,N_AP[:,np.newaxis],T_total)

            # integrate rho_nu with poisson
            logl = np.zeros(params.shape[0])

            #if (loop):
            #    p_N_AP_arr = np.zeros((self.N_AP.shape[0],params.shape[0]))

            binom = sp_binom(self.n,self.k_AP)

            for i,(gamma,delta,nu_max) in enumerate(params):
                #$print('params:',g,d,n)
                fail = False

                if loop:
                    p_N_AP_arr = np.zeros(self.N_AP.shape[0])
                    for j,N in enumerate(self.N_AP):
                        p_N_AP_arr[j],_ = quad(f,0,nu_max,args=(gamma,delta,nu_max,N,self.T_total))
                else:
                    eps_pow = -8
                    while True:
                        if eps_pow==-1:
                            #print('tolerance too high - breaking!')
                            #output = np.ones(inputs[0])*(-100)
                            fail = True
                            logl[i] = -10000
                            #return -10**10 * np.ones(params.shape[0])
                            break
                        try:
                            p_N_AP_arr,err = quadpy.quad(f,[0],[nu_max],
                                args=(gamma,delta,nu_max,self.N_AP,self.T_total),
                                epsabs=10**eps_pow, epsrel=10**eps_pow,
                                limit=1000)
                            break
                        except:
                            eps_pow += 1

                # calculate chance to appear k times
                if (not fail):

                    if loop:
                        p_k = binom * p_N_AP_arr**self.k_AP * (1-p_N_AP_arr)**(self.n - self.k_AP)
                    else:
                        p_k = binom * p_N_AP_arr[:,0]**self.k_AP * (1-p_N_AP_arr[:,0])**(self.n - self.k_AP)

                    #print(p_N_AP_arr.shape,p_k.shape)

                    #print(f"zero counts: {(p_k>0).sum()}")

                    logl[i] = np.log(p_k[p_k>0]).sum()
                    #logl[i] = np.log(p_k[p_k>0]).sum()

            #return np.log(p_k[p_k>0]).sum()
            #print(f"log likelihood:",logl)
            return logl
        return loglikelihood

    ## currently issues: ~99% of traces diverge -> pole is difficult to fit
    def run_on_data(self,draws=5000,tune=10000,loadPath=None,savePath=None,**kwargs):

        """
            ToDo:
                * adjust logp method, such that it calculates
                    \int_0^1/T logp(nu) dnu
                to estimate probability of 0Hz measurements, instead of point-measure at 1/T
                -> use scipy.integrate.quad
                * write documentation!
        """

        if loadPath:
            self.trace = az.from_netcdf(loadPath)
            return

        logl = self.setLogLikelihood(loop=False)
        priorTrafo = self.setPriorTransform()

        for i in range(0,I.data_df.shape[1]):
            time_start = time.time()
            self.prepare_data(i,withZeros=False)
            sampler = ultranest.ReactiveNestedSampler(
                self.param_names,
                logl,
                priorTrafo,
                log_dir=f"logs_animal_{i}",
                vectorized=True
            )
            result = sampler.run(
                show_status=True,
                min_num_live_points=400,
                max_num_improvement_loops=3
            )
            cornerplot(result)
            dt = time.time() - time_start
            print(f"time taken to analyze animal {i}: {dt}s")

        return result, sampler



I = Inference()
I.load_data('empirical',filePath='../../data/BuscheLab/2P_data.xlsx',include_silent=True)
result, sampler = I.run_on_data()
