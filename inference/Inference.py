import os
import numpy as np
import pymc3 as pm
import theano.tensor as tt
import arviz as az

from empirical import *

class Inference:

    dataLoaded = False
    paras = {}

    def __init__(self):

        os.environ['MKL_NUM_THREADS'] = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = '1'

    def set_model(self, func, paras=None):
        assert self.dataLoaded, "Data is not yet loaded, but required for setting the model parameters! (using 'load_data')"

        if callable(func):
            logp = func
        elif func=='selfcon':

            if not paras:
                paras = {
                    'gamma':{'mu':1.5, 'sigma':1.0, 'sigma_animal':1.0,'prior':'Normal'},
                    'delta':{'mu':4.,'sigma':2., 'sigma_animal':1.0,'prior':'Normal'},
                    'nu_max':{'mu':60.,'sigma':20., 'sigma_animal':5.0,'prior':'Normal'}
                }
            assert all(key in paras for key in ("gamma","delta","nu_max")), "Please provide all necessary parameters!"
            def logp(data,paras):
                scaled_NU = tt.log(data / paras['nu_max'])
                return - tt.log( paras['nu_max'] / paras['gamma'] * tt.sqrt( -np.pi * scaled_NU ) ) - paras['delta']**2 / 2 + \
                    ( paras['gamma']**2 - 1 ) * scaled_NU + \
                    tt.log( tt.cosh( paras['gamma'] * paras['delta'] * tt.sqrt( -2 * scaled_NU ) ) )
            # def logp_raw(data,paras):
            #     scaled_NU = np.log(data / paras['nu_max'])
            #     return - np.log( paras['nu_max'] / paras['gamma'] * np.sqrt( -np.pi * scaled_NU ) ) - paras['delta']**2 / 2 + \
            #         ( paras['gamma']**2 - 1 ) * scaled_NU + \
            #         np.log( np.cosh( paras['gamma'] * paras['delta'] * np.sqrt( -2 * scaled_NU ) ) )

        elif func=='lognorm':
            assert all(key in paras for key in ("mu","sigma")), "Please provide all necessary parameters!"

            def logp(data,paras):
                return 0


        self.logp = logp
        # self.logp_raw = logp_raw

        for key in paras.keys():
            self.paras[key] = {}
            for para in paras[key].keys():
                self.paras[key][para] = paras[key][para]


    def load_data(self, dataType='empirical', filePath='../data/BuscheLab/2P_data.xlsx'):

        self.mP = ModelParams(dataType, filePath=filePath, population_keys=['*mouse_type','animal'])
        self.data_df = self.mP.regularize_rates()

        self.data = self.data_df.to_numpy()

        N_zeros = (self.data==0).sum()
        print(f'zeros in data: {N_zeros}')
#        self.data[self.data==0] = np.random.rand(N_zeros)*1./600

        T = 600.
        self.data[self.data<=1./T] = -np.log(1-np.random.rand((self.data<=1./T).sum()))/T

        self.data_mask = ~np.isnan(self.data) & (self.data>0)

        self.dataLoaded = True


    def construct_model_hierarchical(self, name):

        """
            function to create a hierarchical model of the parameter 'name'.

            Input:
                name: string
                    name of the parameter, used for distribution name
                mP: ModelParams Object
                    including 'rates' DataFrame, arranged in a multicolumnar structure
                    reflecting the different populations. Column names starting with '*'
                    are used on the hierarchical level to infer mean value lower level
                    distributions
                mu_population: float
                    mean value of top-level prior
                sigma_population: float
                    std of top-level prior
                sigma_animal:
                    std of low-level prior

            Output:
                prior_animal: pm.distribution
                    distribution containing a draw for each neuron in the dataset
                    ! this should be masked, when populations contain different neuron
                    numbers

            ToDo:
                * enable using different distributions, other than 'Normal' (default)
                * in both, top and bottom level?
                * what about Student-T?
        """

        # obtain information from data-column-header
        population_names = list(self.mP.rates.columns.names)
        N_levels = len(population_names)

        # only columns with '*' are used for top-level inference
        is_hierarchical = [name.startswith('*') for name in population_names]

        print(f'name: {name}')
        print(self.paras[name])

        hierarchical_shape = [
            sha if is_hierarchical[p] else 1
            for p,sha in enumerate(self.mP.data_shape)
        ]
        # print(hierarchical_shape)

        # create top-level distribution for prior mean-values
        prior_mean = self.paras[name]['mu'] + \
            pm.Normal(f'{name}_population',
                mu=0.,
                sigma=self.paras[name]['sigma'],
                shape=hierarchical_shape
            )
        # tt.printing.Print('prior mean shape')(tt.shape(prior_mean))

        # create real distribution, from which values are drawn
        prior_animal = pm.Normal(f'{name}',
            mu=prior_mean,
            sigma=self.paras[name]['sigma_animal'],
            shape=self.mP.data_shape
        )

        # tt.printing.Print('prior animal')(tt.shape(prior_animal))
        prior_animal = prior_animal.reshape((1,-1))

        # tt.printing.Print('prior animal')(tt.shape(prior_animal))

        # broadcast distribution to draw values for each neuron
        prior_animal = tt.tile(
            prior_animal,
            (self.data.shape[0], 1) ### why -1?    *[1]*(N_levels-1)
        )
        tt.printing.Print('prior animal (final)')(tt.shape(prior_animal))

        return prior_animal


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
            return az.from_netcdf(loadPath)

        data_observed = self.data[self.data_mask]
        #silent = tt.le(data_observed,1./600)
        with pm.Model() as model:
            # replace normals with student-t distributions

            priors = {}
            for para in self.paras:
                priors[para] = self.construct_model_hierarchical(para)
                priors[para] = priors[para][self.data_mask]


            def likelihood(data):

                # introduce checks for consistency, etc
                logP = self.logp(data,priors)

                # penalize nan-entries (e.g. when log is negative, etc)
                # logP_masked = tt.switch(tt.isnan(logP), 0, logP)
                # logP = tt.switch(silent, logP-10., logP)
                #min_val = tt.min(logP_masked)
                #tt.printing.Print('logP minimum')(tt.min(logP_masked))
                tt.printing.Print('logP')(logP)

                #logP = tt.switch(tt.isnan(logP), min_val*2, logP)

                return tt.sum(logP)

            ## watch out: for some reason, NaNs in observed data are converted to 0s
            logP = pm.DensityDist('logP',likelihood,observed=data_observed)

            trace = pm.sample(
                init='adapt_diag',
                step=pm.Metropolis(),
                chains=4,draws=draws,tune=tune,
                return_inferencedata=True,
                **kwargs)

            if savePath:
                trace.to_netcdf(savePath)

            self.trace = trace
