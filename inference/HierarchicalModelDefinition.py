import numpy as np
import inspect
import pprint

class HierarchicalModel:
    """
        Defines a general class for setting up a hierarchical model for bayesian inference. Has to be inherited by a specific model class, which then further specifies the loglikelihood etc
    """

    def __init__(self,nSamples):

        self.nSamples = nSamples

    def set_priors(self, priors_init):
        """
            Set the priors for the model. The priors are defined in the priors_init dictionary, which has to follow the following structure:

            priors_init = {
                **param1** : {
                    'hierarchical': {
                        'function':     **transform_function**,
                        'params': {
                            'loc':      'mean',
                            'scale':    **ref_to_2nd_key**
                        },
                    },
                    mean : {
                        'function':     **transform_function**,
                        'params':       {**params_for_function**}
                    },
                    **2nd_key** : {
                        'function':     **transform_function**,
                        'params':       {**params_for_function**}
                    }
                },
                **param2** : ...,
            }
            with **x** being placeholders for the actual values.

            All parameters appearing in "hierarchical" will be treated as hierarchical parameters. If a parameter is hierarchical, the mean and **2nd_key** parameter define the meta distribution. If a parameter is not hierarchical, only the mean parameter is used.
        """

        self.paramNames = []
        self.paramIn = list(priors_init.keys())
        self.priors = {}

        ct = 0
        for prior_key, prior in priors_init.items():

            if prior.get("hierarchical"):

                ## add the parameters for the hierarchical prior
                for sub_key, sub_prior in prior["parameters"].items():

                    self.set_prior_param(
                        sub_prior,
                        ct,
                        prior_key,
                        sub_key,
                    )
                    ct += 1

                ## then, add the actual parameters for the hierarchical prior
                self.set_prior_param(prior, ct, prior_key, lower_hierarchy_level=True)
                ct += self.nSamples
            # if prior["function"] is None:

            else:
                ## add the parameters for the non-hierarchical prior
                self.set_prior_param(prior, ct, prior_key)
                ct += 1 if prior["sample"] else 0

        self.nParams = len(self.paramNames)
        self.wrap = np.zeros(self.nParams).astype('bool')

    def set_prior_param(
        self, priors_init, ct, param, key=None, lower_hierarchy_level=False
    ):
        """
        TODO
        * description of the function
        """

        paramName = param + (f"_{key}" if key else "")
        self.priors[paramName] = {}

        self.priors[paramName]["n"] = self.nSamples if lower_hierarchy_level else 1
        self.priors[paramName]["idx"] = ct

        if priors_init["function"] is None:
            # pprint.pprint(priors_init)
            self.priors[paramName]["value"] = list(priors_init["parameters"].values())[
                0
            ]
            self.priors[paramName]["transform"] = None

        elif lower_hierarchy_level:
            for i in range(self.nSamples):
                self.paramNames.append(f'{param}_{i}')

            # get indexes of hierarchical parameters for quick access later
            self.priors[paramName]["paramNames"] = []
            for var in priors_init["parameters"].keys():
                self.priors[paramName][f"idx_{var}"] = self.priors[f"{param}_{var}"][
                    "idx"
                ]
                self.priors[paramName]["paramNames"].append(var)

            self.priors[paramName]["transform"] = lambda x, params, fun=priors_init[
                "function"
            ]: fun(x, **params)

            # self.priors[paramName]["paramNames"] = priors_init["parameters"].keys()
        else:
            self.paramNames.append(paramName)

            self.priors[paramName]["transform"] = (
                lambda x, params=priors_init["parameters"], fun=priors_init[
                    "function"
                ]: fun(x, **params)
            )

    def set_prior_transform(self,vectorized=True):
        '''
            sets the prior transform function for the model

            only takes as input the mode, which can be either of
            - 'vectorized': vectorized prior transform function
            - 'scalar': scalar prior transform function
            - 'tensor': tensor prior transform function
        '''

        def prior_transform(p_in):

            """
                The actual prior transform function, which transforms the random variables from the unit hypercube to the actual priors
            """

            if len(p_in.shape)==1:
                p_in = p_in[np.newaxis,...]
            p_out = np.zeros_like(p_in)

            for prior in self.priors.values():
                # print(key, prior)
                if prior["transform"] is None:
                    continue

                if prior["n"] == 1:

                    p_out[:, prior["idx"]] = prior["transform"](p_in[:, prior["idx"]])

                else:
                    params = {}

                    for var in prior["paramNames"]:
                        params[var] = p_out[:, prior[f"idx_{var}"], np.newaxis]

                    p_out[:, prior["idx"] : prior["idx"] + prior["n"]] = prior[
                        "transform"
                    ](p_in[:, prior["idx"] : prior["idx"] + prior["n"]], params=params)
            # print('proposed:',p_out)
            if vectorized:
                return p_out
            else:
                return p_out[0,:]

        return prior_transform

    def get_params_from_p(self, p_in, idx_chain=None, idx_sample=None):
        """
        obtains a readable dictionary of parameters from the input p_in
        """
        params = {}
        nChains = p_in.shape[0]

        assert (
            idx_chain is None or idx_chain < nChains
        ), "idx_chain must be smaller than the number of chains in p_in"
        assert (
            idx_sample is None or idx_sample < self.nSamples
        ), "idx_sample must be smaller than the number of samples in p_in"

        slice_chain = slice(None) if idx_chain is None else idx_chain

        for var in self.paramIn:

            params[var] = np.zeros(
                ((nChains,) if idx_chain is None else ())
                + ((self.nSamples,) if idx_sample is None else ())
            )

            if self.priors[var].get("transform"):

                if idx_sample is None:
                    slice_sample = slice(
                        self.priors[var]["idx"],
                        self.priors[var]["idx"]
                        + (self.nSamples if self.priors[var]["n"] > 1 else 1),
                    )
                else:
                    p_idx_sample = self.priors[var]["idx"] + (
                        idx_sample if self.priors[var]["n"] > 1 else 0
                    )
                    slice_sample = slice(p_idx_sample, p_idx_sample + 1)

                # Get the sliced values from p_in
                sliced = np.squeeze(p_in[slice_chain, slice_sample])
            else:
                sliced = np.squeeze(self.priors[var]["value"][slice_chain])

            # Fill params[var] with the sliced values, handling both scalar and array cases
            if params[var].shape == ():
                params[var] = sliced
            else:
                params[var][...] = sliced

        return params
        # if idx_chain is None and idx_sample is None:
        # if no specific chain or sample is requested, return all
        # else:
        # return np.squeeze(params)


def prior_structure(function=None, **kwargs):
    """
    creates a dictionary with appropriate structure for the prior distribution

    inputs:
        function:   None | callable
            None:       no sampling performed, only value in params is returned
            callable:   function to be applied to the parameters
        kwargs:     dict
            key-value pairs of parameters to be passed to the function
    """

    if function is None:
        ### watch out: name of variable is arbitrary!
        if len(kwargs) != 1:
            raise ValueError(
                "If function is None, kwargs must contain exactly one entry."
            )
    else:
        if not callable(function):
            raise ValueError("function must be callable if not None.")
        # Try to check if kwargs match the function signature
        sig = inspect.signature(function)
        params = set(sig.parameters.keys())
        missing = params - set(kwargs.keys()) - {"x"}
        if missing:
            raise ValueError(f"Missing required kwargs for function: {missing}")

    return {
        "hierarchical": (
            np.any([isinstance(val, dict) for val in kwargs.values()])
            if kwargs
            else False
        ),
        "sample": not (function is None),
        "function": function,
        ## could also just contain the params as entries here?!
        # **kwargs,
        "parameters": kwargs,
    }
