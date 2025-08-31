import logging
import numpy as np
import re

# from typing import List, Tuple
# import pprint
from matplotlib import pyplot as plt


from scipy.special import gammaln

from .HierarchicalBayesModel import HierarchicalModel, prior_structure, norm_ppf, halfnorm_ppf

from collections import Counter


from .utils.utils import adaptive_integration, spike_observation_probability

# from inference.transform_meta_to_bio import get_nu_bar, get_q, get_tau_I, get_alpha_0
# from inference.transform_bio_to_meta import get_nu_max, get_gamma, get_delta

from inference.network import Network

from .utils.structures import DistributionModelParams as distr_params
from .utils.utils import rho_nu

import warnings
warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO)


class BayesModel(HierarchicalModel):
    
    def prepare_data(self, event_counts, T, **kwargs):

        """
            Prepares the data for the model by setting the relevant dimensions 
            for shape and obtaining the number of populations from input shape
        """
        dims = event_counts.shape
        iter_dims = np.ones_like(dims, dtype=bool)
        iter_dims[-2:] = False


        super().prepare_data(event_counts, T, iter_dims=iter_dims, **kwargs)
        self.dimensions["n_pop"] = dims[-2]   ## obtain number of populations

    def set_logl(
        self,
        vectorized=False,
        correct_N=0,
        bias_to_expected_max=0,
        bias_to_mean=0,
        correct_threshold=10 ** (-4),
        biological=False,
    ):
        """
        sets the log likelihood function for the model,
        which calculates the log likelihood of the data given the model parameters

        it assumes that neuron activity is provided as spike counts per neuron
        and the model is a Poisson model

        input
            vectorized [bool]
                whether the log likelihood function should be vectorized or not
            correct_N [int]
                spike number up to which the likelihood should be corrected, if required

        """

        offset = 0.5

        self.calculate_binom()
        
        print("Different values of kappa currently not possible - check what needs to be done!")

        if biological:
        # initialize the network class with according populations & synapses
            self.net = Network()
            for key in self.parameter_names:

                var,(p,s) = parse_name_and_indices(key, ["pop","s"])
                if p is not None:
                    self.net.register_population(p,J0=(-1.)**(p+1))

                if p is not None and s is not None:
                    self.net.populations[p].register_synapse(s)



        def loglikelihood(p_in):
            """
            log likelihood function for the model

            The likelihood is calculated as count based for low expected spike counts (nu_max*T < N_low_counts)
            and as continuous for medium and high expected spike counts (nu_max*T >= N_low_counts)

            input (* = class-variable)
                p_in [nChain, nAnimals, nParams]
                    The input parameters for the model, which are the gamma, delta and nu_max

            output
                logp [nChain]
                    The log likelihood of the data, given the model parameters
                        sum_{nAnimals} sum_{N_max} log( p_k_AP )

            """
            # define dimensions
            if len(p_in.shape) == 1:
                p_in = p_in[np.newaxis, :]

            ## should be done somewhere else
            # max_spike_count = np.nanmax(self.event_counts, axis=-1).astype("int")

            n_chain = p_in.shape[0]
            logl = np.zeros((n_chain,) + self.dimensions["shape_iter"])

            for idx in self.dimensions["iterator"]:
                
                # cts = self.data["event_counts"][idx]
                max_spike_count = int(np.nanmax(self.data["event_counts"][idx]))
                
                for chain in range(n_chain):
                    """
                        get the parameters from the input p_in into dictionary
                    """

                    full_idx = (chain,) + idx

                    params = self.get_params_from_p(p_in, idx_chain=chain, idx=idx, biological=biological)
                    if not params:
                    # happens if parameters are weird
                        logl[full_idx] = -(10**6)
                        continue

                    ## calculate maximum number of spikes:
                    max_spike_count_model = np.squeeze(
                        params["distr"][0].nu_max * self.T
                    ).astype("int")
                    # print(f"{max_spike_count_model=} vs {max_spike_count=}")

                    logl[full_idx], abort = self.parameter_penalties(params, max_spike_count_model, idx)
                    if abort: 
                        continue

                    """
                        from here, animal dependent calculations necessary
                        (due to differences in nu_max and thus max_spike_count_model)
                    """
                    for p in range(self.dimensions["n_pop"]):
                        idx_population = idx + (p,)
                        ## remove all spike counts that are higher than max_spike_count_model
                        N_AP_empirical = self.binom[idx_population]["N_AP"][
                            self.binom[idx_population]["N_AP"] < max_spike_count_model
                        ]
                        k_AP_empirical = self.binom[idx_population]["k_AP"][
                            self.binom[idx_population]["N_AP"] < max_spike_count_model
                        ]

                        N_AP = np.arange(max_spike_count_model).astype("int")
                        k_AP = np.zeros(max_spike_count_model, dtype="int")
                        k_AP[N_AP_empirical] = k_AP_empirical

                        ### calculate actual log-likelihood
                    # print(params)
                        p_N_AP = get_p_N_AP(
                            (N_AP + offset) / self.T,
                            (params["distr"][p],),
                            self.T,
                            correct_N=correct_N,
                            correct_threshold=correct_threshold,
                        )
                        # print("p_N_AP", p_N_AP[:5])

                        # plt.figure()
                        # plt.plot(N_AP, p_N_AP, label=f"Population {p}")

                        log_binom = self.binom[idx_population]["log_binomial_coefficient"][k_AP]

                        # logl_binom = (
                        #     log_binom
                        #     + k_AP * np.log(p_N_AP[N_AP])
                        #     + (self.data["n_neurons"][idx_population] - k_AP)
                        #     * np.log(1 - p_N_AP[N_AP])
                        # )
                        # plt.plot(N_AP, logl_binom, label=f"Binom Pop {p}")
                        # plt.xlabel("Number of Action Potentials")
                        # plt.ylabel("Probability")
                        # plt.title("Probability Distribution of Action Potentials")
                        # plt.legend()
                        # plt.show()


                        # print("logl (pre binom)",logl[full_idx])

                        logl[full_idx] += (
                            log_binom
                            + k_AP * np.log(p_N_AP[N_AP])
                            + (self.data["n_neurons"][idx_population] - k_AP)
                            * np.log(1 - p_N_AP[N_AP])
                        ).sum()

                        # print("logl (post binom)",logl[full_idx])


                        ### add penalties
                        # print("penalties:", self.penalty_nu_max(max_spike_count_model, animal))

                        # if bias_to_mean > 0:
                        #     logl[full_idx] += self.penalty_mean_bias(
                        #         params, bias_to_mean, idx
                        #     )


                        if bias_to_expected_max > 0 and (
                            max_spike_count < max_spike_count_model
                        ):

                            p_max = expected_maximum_value(
                                p_N_AP,
                                self.data["n_neurons"][idx],
                                max_spike_count,
                            )

                            logl[full_idx] += np.log(p_max) * bias_to_expected_max

                        elif bias_to_expected_max > 0:
                            logl[full_idx] += bias_to_expected_max * -(10**6)

                    # self.log.debug((f'logls: {logl_low=}, {logl_intermediate=}, {logl_tail=}'))
                    if not np.isfinite(logl[full_idx]):
                        logl[full_idx] = -(10**6)

            if vectorized:
                return logl.sum(axis=tuple(range(1, self.dimensions["n"])))
            else:
                self.log.debug(("logl:", logl[0, :].sum()))
                return logl[0, :].sum()

        return loglikelihood

    def calculate_binom(self):

        ## get the maximum number of spikes
        max_spike_count = np.nanmax(self.data["event_counts"]).astype("int64")

        ## calculate the number of occurences 'k_AP' of each spike number per animal,
        ## as well as the binomial coefficient of occurence 'k_AP', given n neurons
        self.binom = {}

        for idx in self.dimensions["iterator"]:

            for p in range(self.dimensions["n_pop"]):
                idx_p = idx + (p,)
                event_counts = self.data["event_counts"][idx_p]
                # print("count shape:",event_counts.shape)
                N = event_counts[np.isfinite(event_counts)].astype("int64")
                N_ct = Counter(N)

                k_AP_all = np.zeros(max_spike_count + 1)
                k_AP_all[list(N_ct.keys())] = list(N_ct.values())

                # Store results in binom dict
                self.binom[idx_p] = {}
                self.binom[idx_p]["N_AP"] = np.where(k_AP_all > 0)[0]
                self.binom[idx_p]["k_AP"] = k_AP_all[self.binom[idx_p]["N_AP"]].astype("int")

                self.binom[idx_p]["log_binomial_coefficient"] = log_binomial(
                    self.data["n_neurons"][idx_p],
                    np.arange(self.data["n_neurons"][idx_p] + 1),
                )


    def parameter_penalties(self,params, max_spike_count_model, idx):
                        
        # for pop in params["distr"]:
        #     if not pop.are_values_ok():
        #         return -(10**6), True

        ## calculate penalty for model max spike count less than observed data
        logl_penalty = 0

        for p in range(self.dimensions["n_pop"]):
            idx_population = idx + (p,)
            idx_too_high = np.where(self.binom[idx_population]["N_AP"] > max_spike_count_model)[0]
            for i in idx_too_high:
                ## get all spike counts higher than allowed by model
                N_AP_too_high = self.binom[idx_population]["N_AP"][i]
                k_AP_too_high = self.binom[idx_population]["k_AP"][i]
                logl_penalty += (
                    -(10**6)
                    * k_AP_too_high
                    * ((N_AP_too_high - max_spike_count_model) / max_spike_count_model) ** 2
                )
        # print(f"logl_penalty: {logl_penalty}")


        return logl_penalty, False


    # def penalty_mean_bias(self, params, bias_to_mean, idx):

    #     if self.two_pop:
    #         nu_mean = params["p"] * get_nu_bar(**params["distr"][0])
    #         nu_mean += (1 - params["p"]) * get_nu_bar(**params["distr"][1])
            
    #         nu_SD = params["p"] * get_q(**params["distr"][0])
    #         nu_SD += (1 - params["p"]) * get_q(**params["distr"][1])
    #     else:
    #         nu_mean = get_nu_bar(**params["distr"][0])
    #         nu_SD = get_q(**params["distr"][0])

    #     nu_sigma = np.sqrt(nu_SD - nu_mean**2) / np.sqrt(self.data["n_neurons"][idx])

    #     rates = (
    #         self.data["event_counts"][
    #             idx, np.isfinite(self.data["event_counts"][idx, :])
    #         ]
    #         / self.T
    #     )
    #     logl_bias_to_mean = -((rates.mean() - nu_mean) ** 2) / (2 * nu_sigma**2)

    #     return logl_bias_to_mean * bias_to_mean
    #     # logl[i, a] += bias_to_mean * logl_bias_to_mean


    def get_params_from_p(self, p_in, idx_chain=None, idx=None, biological=False):
        
        """
            build structure of distribution parameters from p_in to dictionary shape as required by logl
        """

        params_tmp = super().get_params_from_p(p_in, idx_chain=idx_chain, idx=idx)

        if biological:

            for key,val in params_tmp.items():

                var,(p,s) = parse_name_and_indices(key,["pop","s"])

                print("set params:")
                print(key,val)
                print(var,p,s)

                if p is None and s is None:
                    setattr(self.net,var,val)
                elif p is not None and s is None:
                    setattr(self.net.populations[p],var,val)
                elif p is not None and s is not None:
                    setattr(self.net.populations[p].synapses[s],var,val)
                else:
                    print(f"Unexpected parameter structure for {key}: {val}")
            
            is_ok = self.net.are_values_ok()
            if not is_ok:
                return False

            self.net.set_weights()
            self.net.calculate_sigma_V()
            self.net.solve_selfcon()

            params = {}
            for p,pop in enumerate(self.net.populations):
                for key in ["gamma","delta","nu_max"]:
                    params[f"{key}_pop{p}"] = getattr(pop, key)

            # print("flat:",params)
            params = build_distr_structure_from_params(params, self.parameter_names)
            # print("full:",params)

        else:

            params = build_distr_structure_from_params(params_tmp, self.parameter_names)

        # if len(params["distr"])==2: #self.two_pop:
        #     params["distr"][1]["nu_max"] = params["distr"][0]["nu_max"]

        self.log.debug("\n params:", params)

        return params


def build_distr_structure_from_params(params_tmp, paramNames):

    # print(params_tmp)
    paramNames = params_tmp.keys()#["gamma", "delta", "nu_max"]
    params = {"distr": []}
    for var in paramNames:
        var_split = var.split("_")
        var_root = "_".join(var_split[:-1])
        idx_population = int(var_split[-1][-1]) if var_split[-1].startswith("pop") else np.nan
        if len(params["distr"]) < idx_population + 1:
            params["distr"].append(distr_params(0,0,0))

        # print(var,var_root,idx_population)
        if np.isnan(idx_population):
            # print("when does this show up?")
            # print(f"{var=}, {var_root=}, {idx_population=}")
            params[var] = params_tmp[var]
        else:
            setattr(params["distr"][idx_population],var_root,params_tmp[var])

    return params


def expected_maximum_value(p, n, max_value=None):
    """
        obtain distribution of extreme values, given an input distribution p and a sample size n.
        If max_value is provided, return the probability of observing that value

    p   array(float)
        - input distribution (probability mass function)
    n   int
        - sample size

    output:
        - expected maximum value probability distribution
    """

    ### now, get extreme value distribution
    p_cum = np.pad(
        np.cumsum(p) ** n,
        (1, 0),
        mode="constant",
        constant_values=0,
    )
    p_extreme = np.diff(p_cum)
    p_extreme /= p_extreme.sum()  # normalizing

    if max_value is None:
        return p_extreme
    else:
        return p_extreme[max_value]



def get_p_N_AP(nu, args_rho, T, correct_N=5, correct_threshold=10 ** (-4)):

    # print(f"{args_rho=}")
    p_N_AP = rho_nu(nu, *args_rho) / T
    p_N_AP[~np.isfinite(p_N_AP)] = np.nan

    if correct_N == 0:
        return p_N_AP
    
    ## correcting terms until normalization becomes decent
    p_full = p_N_AP[:correct_N].copy()

    for N in range(correct_N):
        ### correcting one at a time
        p_N_AP[N] = adaptive_integration(spike_observation_probability, 0, 100.0 / T, args=(args_rho,(N, T)), eps_pow=-2)

        idx_start = max(0, N - 3)
        dp = np.abs(p_full[idx_start : N + 1] - p_N_AP[idx_start : N + 1])

        ## ensure that p_N_AP converges and is not just oscillating
        if np.all(dp < correct_threshold):  # all recent corrections have been small
            # print(f"approximation good enough @{N=}, {dp[-1]}")
            break

    return p_N_AP


def get_default_priors():

    prior = {}
    prior["gamma_0"] = prior_structure(norm_ppf, mean=2.0, sigma=0.5)
    prior["delta_0"] = prior_structure(norm_ppf, mean=6.0, sigma=2.0)
    prior["nu_max_0"] = prior_structure(
        None,
        value=25.0,
    )
    return prior

"""
    could be moved to other file?!
"""

import re
from typing import Iterable, List, Optional, Tuple

def parse_name_and_indices(s: str, literals: Iterable[str]) -> Tuple[str, List[Optional[int]]]:
    """
    Returns (variable_name, [idx_or_None per literal in the same order]).
    Variable name = prefix before the first <literal><digits> token.
    """
    lits = list(literals)
    alts = "|".join(map(re.escape, lits))
    # Don't match inside letter-words; allow underscores and punctuation as separators.
    rx = re.compile(rf"(?<![A-Za-z])({alts})(\d+)(?![A-Za-z])")

    found = {}
    first_pos = None

    for m in rx.finditer(s):
        lit, num = m.group(1), int(m.group(2))
        if lit not in found:  # keep only first per literal
            found[lit] = (m.start(), num)
            if first_pos is None or m.start() < first_pos:
                first_pos = m.start()

    name = s[:first_pos-1] if first_pos is not None else s
    indices: List[Optional[int]] = [found[lit][1] if lit in found else None for lit in lits]
    return name, indices



# def get_idx_from_key(key, prefixes):
#     # Search for patterns "prefixX" in key and return X
#     for prefix in prefixes:
#         match = re.search(rf'{prefix}(\d+)', key)
#     return int(match.group(1)) if match else None


"""
    could be moved to helper file
"""
def log_binomial(n, k):
    return gammaln(n + 1) - gammaln(n - k + 1) - gammaln(k + 1)


# def weight(x, mode="sigmoid", **kwargs):
#     # return a + (1-a)/np.log(np.exp(1)+1/(1-x))
#     if mode == "sigmoid":
#         return kwargs["offset"] + (1 - kwargs["offset"]) / (
#             1 + np.exp(kwargs["slope"] * (x - kwargs["threshold"]))
#         )
#     else:
#         return None


# # def weight(x, a=0.5, b=10):
# #     return a + (1 - a) * (1 - x) ** b