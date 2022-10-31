import pymc3 as pm
import numpy as np
import os
import theano.tensor as tt
import arviz as az
from .models import *

## currently issues: ~99% of traces diverge -> pole is difficult to fit
def run_on_data(mP,draws=5000,tune=10000,loadPath=None,savePath=None):

    if loadPath:
        return az.from_netcdf(loadPath)

    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    with pm.Model() as model:
        # replace normals with student-t distributions

        gamma = construct_model_hierarchical('gamma',mP);
        delta = construct_model_hierarchical('delta',mP);
        nu_max = construct_model_nu_max('nu_max',mP);

        def likelihood(spike_data):

            # introduce checks for consistency, etc
            scaled_NU = tt.log(spike_data / nu_max)
            logP = - tt.log( nu_max / gamma * tt.sqrt( -np.pi * scaled_NU ) ) - delta**2 / 2 + \
                ( gamma**2 - 1 ) * scaled_NU + \
                tt.log( tt.cosh( gamma * delta * tt.sqrt( -2 * scaled_NU ) ) )

            # penalize nan-entries (e.g. when log is negative, etc)
            logP_masked = tt.switch(tt.isnan(logP), 0, logP)
            min_val = tt.min(logP_masked)
            tt.printing.Print('logP minimum')(tt.min(logP_masked))

            logP = tt.switch(tt.isnan(logP), min_val*2, logP)

            return tt.sum(logP)

        ## watch out: for some reason, NaNs in observed data are converted to 0s
        logP = pm.DensityDist('logP',likelihood,observed=mP.spikes)

        trace = pm.sample(
            init='adapt_diag',
            chains=4,draws=draws,tune=tune,
            target_accept=0.8,
            return_inferencedata=True)

        if savePath:
            trace.to_netcdf(savePath)

        return trace
