import os, logging
import numpy as np

import pprint
import itertools
from scipy.special import binom as sp_binom, erfinv, gammaln

from .HierarchicalBayesModel import HierarchicalModel, prior_structure

from collections import Counter

from matplotlib import pyplot as plt


from .utils.utils import p_nu, adaptive_integration, f

from DM_theory.functions import get_nu_bar, get_q, get_tau_I, get_alpha_0
from DM_theory.functions import get_nu_max, get_gamma, get_delta

import warnings
warnings.filterwarnings("ignore")


logging.basicConfig(level=logging.INFO)


class BayesModel(HierarchicalModel):

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
            withZeros [bool]
                whether the likelihood should be calculated for spike counts of 0
            correct_N [int]
                spike number up to which the likelihood should be corrected, if required

        """

        offset = 0.5

        self.calculate_binom()

        def loglikelihood(p_in, plot=False):
            """
            log likelihood function for the model

            The likelihood is calculated as count based for low expected spike counts (nu_max*T < N_low_counts)
            and as continuous for medium and high expected spike counts (nu_max*T >= N_low_counts)

            input (* = class-variable)
                p_in [nChain, nAnimals, nParams]
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
            """
            # define dimensions
            if len(p_in.shape) == 1:
                p_in = p_in[np.newaxis, :]

            ## should be done somewhere else
            # max_spike_count = np.nanmax(self.event_counts, axis=-1).astype("int")

            n_chain = p_in.shape[0]

            logl = np.zeros((n_chain,) + self.dimensions["shape"][:-1])

            # print(self.dimensions)

            for idx in self.dimensions["iterator"]:
                # print(idx)

                # for animal in range(self.dims["n_samples"]):

                # print(self.data["event_counts"].shape)
                cts = self.data["event_counts"][idx]
                # print(idx, cts.shape)
                # rates = self.rates[np.isfinite(self.rates[:, animal]), animal]
                max_spike_count = int(np.nanmax(self.data["event_counts"][idx]))
                # print(f"{max_spike_count=}")

                for chain in range(n_chain):

                    """
                    get the parameters from the input p_in into dictionary
                    """

                    full_idx = (chain,) + idx
                    params = self.get_params_from_p(p_in, idx_chain=chain, idx=idx)

                    # print(f"params: {params}")

                    if (
                        np.isnan(params["distr"][0]["nu_max"])
                        or np.isnan(params["distr"][0]["gamma"])
                        or np.isnan(params["distr"][0]["delta"])
                        or params["distr"][0]["nu_max"] < 0
                        or params["distr"][0]["gamma"] < 0
                        or params["distr"][0]["delta"] < 0
                    ):
                        logl[full_idx] = -(10**6)
                        continue

                    if self.two_pop:
                        params["distr"][1]["nu_max"] = params["distr"][0]["nu_max"]

                    ## calculate maximum number of spikes:
                    # print(self.T, params["distr"][0]["nu_max"] * self.T)
                    max_spike_count_model = np.squeeze(
                        params["distr"][0]["nu_max"] * self.T
                    ).astype("int")

                    """
                        from here, animal dependent calculations necessary
                        (due to differences in nu_max and thus max_spike_count_model)
                    """
                    ## remove all spike counts that are higher than max_spike_count_model
                    N_AP_empirical = self.binom[idx]["N_AP"][
                        self.binom[idx]["N_AP"] < max_spike_count_model
                    ]
                    k_AP_empirical = self.binom[idx]["k_AP"][
                        self.binom[idx]["N_AP"] < max_spike_count_model
                    ]

                    N_AP = np.arange(max_spike_count_model).astype("int")
                    k_AP = np.zeros(max_spike_count_model, dtype="int")
                    k_AP[N_AP_empirical] = k_AP_empirical

                    ### calculate actual log-likelihood
                    p_N_AP = get_p_N_AP(
                        (N_AP + offset) / self.T,
                        params,
                        self.T,
                        correct_N=correct_N,
                        correct_threshold=correct_threshold,
                    )

                    log_binom = self.binom[idx]["log_binomial_coefficient"][k_AP]

                    logp = (
                        log_binom
                        + k_AP * np.log(p_N_AP[N_AP])
                        + (self.data["n_neurons"][idx] - k_AP)
                        * np.log(1 - p_N_AP[N_AP])
                    )

                    if plot:
                        fig = plt.figure()
                        ax1 = fig.add_subplot(211)
                        ax1.plot(N_AP, k_AP, label="k_AP")
                        ax2 = fig.add_subplot(212)
                        ax2.axhline(0, color="black", linestyle="--")
                        ax2.plot(N_AP, logp, label="p_N_AP")
                        plt.show(block=False)
                    # * weight(
                    #     N_AP / max_spike_count_model, "sigmoid", offset=0.2, slope=50, threshold=0.1
                    # )

                    logl_measures = logp.sum()
                    logl[full_idx] += logl_measures

                    ### add penalties
                    # print("penalties:", self.penalty_nu_max(max_spike_count_model, animal))
                    logl[full_idx] += self.penalty_nu_max(max_spike_count_model, idx)

                    if bias_to_mean > 0:
                        logl[full_idx] += self.penalty_mean_bias(
                            params, bias_to_mean, idx
                        )

                    if bias_to_expected_max > 0 and (
                        max_spike_count < max_spike_count_model
                    ):
                        # p_extreme = 1 - 1.0 / (params["distr"][0]["nu_max"] - N_AP / self.T + 1) ** 2

                        p_max = expected_maximum_value(
                            p_N_AP,
                            self.data["n_neurons"][idx],
                            max_spike_count,
                        )
                        # p_max *= (
                        #     1.0 - 1.0 / ((max_spike_count_model - max_spike_count) / self.T + 1) ** 2
                        # )

                        logl[full_idx] += np.log(p_max) * bias_to_expected_max

                        # logl[i, animal] -= self.penalty_expected_max_bias(
                        #     params,
                        #     p_N_AP,
                        #     N_AP,
                        #     max_spike_count,
                        #     bias_to_expected_max,
                        #     animal,
                        # )
                    elif bias_to_expected_max > 0:
                        logl[full_idx] += bias_to_expected_max * -(10**6)

                    """
                        sum up the log likelihoods to obtain the overall estimate for this animal
                    """
                    # self.log.debug((f'logls: {logl_low=}, {logl_intermediate=}, {logl_tail=}'))

                    # logl[i,a] = logl_low + logl_intermediate + logl_tail
                    # print(logl[full_idx])
                    if not np.isfinite(logl[full_idx]):
                        logl[full_idx] = -(10**6)

            if vectorized:
                # print(logl.shape, logl)
                # print(tuple(range(1, self.dimensions["n"])))
                return logl.sum(axis=tuple(range(1, self.dimensions["n"])))
            # return logl.sum(axis=1)
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

        # print(self.dimensions)
        # Iterate through all but the last dimension of event_counts
        for idx in self.dimensions["iterator"]:

            event_counts = self.data["event_counts"][idx]
            N = event_counts[np.isfinite(event_counts)].astype("int64")
            N_ct = Counter(N)

            k_AP_all = np.zeros(max_spike_count + 1)
            k_AP_all[list(N_ct.keys())] = list(N_ct.values())

            # Store results in binom dict using idx as key

            self.binom[idx] = {}
            # d = self.binom
            # for i in idx[:-1]:
            #     d = d.setdefault(i, {})
            # print(f"idx: {idx}, d: {self.binom}")

            self.binom[idx] = {}
            self.binom[idx]["N_AP"] = np.where(k_AP_all > 0)[0]
            self.binom[idx]["k_AP"] = k_AP_all[self.binom[idx]["N_AP"]].astype("int")

            self.binom[idx]["log_binomial_coefficient"] = log_binomial(
                self.data["n_neurons"][idx],
                np.arange(self.data["n_neurons"][idx] + 1),
            )
        # print(self.binom)

        # for a, event_counts_animal in enumerate(self.data["event_counts"]):
        #     # print(f"{a=}")
        #     # print(N.shape, N)
        #     self.binom[a] = {}
        #     for c, event_counts_condition in enumerate(event_counts_animal):
        #         N = event_counts_condition[np.isfinite(event_counts_condition)].astype(
        #             "int64"
        #         )
        #         N_ct = Counter(N)

        #         k_AP_all = np.zeros(max_spike_count + 1)
        #         k_AP_all[list(N_ct.keys())] = list(N_ct.values())

        #         self.binom[a][c] = {}

        #         self.binom[a][c]["N_AP"] = np.where(k_AP_all > 0)[0]
        #         self.binom[a][c]["k_AP"] = k_AP_all[self.binom[a][c]["N_AP"]].astype(
        #             "int"
        #         )
        #         self.binom[a][c]["log_binomial_coefficient"] = log_binomial(
        #             self.dims["n_neurons"][a, c],
        #             np.arange(self.dims["n_neurons"][a, c] + 1),
        #         )
        # print(self.binom)

    def penalty_nu_max(self, max_spike_count_model, idx):

        ## get all spike counts higher than allowed by model
        idx_too_high = np.where(self.binom[idx]["N_AP"] > max_spike_count_model)[0]
        logp_penalty = 0
        for idx in idx_too_high:
            # print(idx)
            N_AP_too_high = self.binom[idx]["N_AP"][idx]
            k_AP_too_high = self.binom[idx]["k_AP"][idx]
            # print(f'{N_AP_too_high=}, {k_AP_too_high=}')
            logp_penalty += (
                -(10**6)
                * k_AP_too_high
                * ((N_AP_too_high - max_spike_count_model) / max_spike_count_model) ** 2
            )
            # logp_penalty += -(10**0) * k_AP_too_high * (N_AP_too_high - max_spike_count_model) ** 2
        return logp_penalty

    def penalty_mean_bias(self, params, bias_to_mean, idx):

        if self.two_pop:
            nu_mean = params["p"] * get_nu_bar(
                params["distr"][0]["gamma"],
                params["distr"][0]["delta"],
                params["distr"][0]["nu_max"],
            )
            nu_mean += (1 - params["p"]) * get_nu_bar(
                params["distr"][1]["gamma"],
                params["distr"][1]["delta"],
                params["distr"][1]["nu_max"],
            )
            nu_SD = params["p"] * get_q(
                params["distr"][0]["gamma"],
                params["distr"][0]["delta"],
                params["distr"][0]["nu_max"],
            )
            nu_SD += (1 - params["p"]) * get_q(
                params["distr"][1]["gamma"],
                params["distr"][1]["delta"],
                params["distr"][1]["nu_max"],
            )
        else:
            nu_mean = get_nu_bar(
                params["distr"][0]["gamma"],
                params["distr"][0]["delta"],
                params["distr"][0]["nu_max"],
            )
            nu_SD = get_q(
                params["distr"][0]["gamma"],
                params["distr"][0]["delta"],
                params["distr"][0]["nu_max"],
            )

        nu_sigma = np.sqrt(nu_SD - nu_mean**2) / np.sqrt(self.data["n_neurons"][idx])

        rates = (
            self.data["event_counts"][
                idx, np.isfinite(self.data["event_counts"][idx, :])
            ]
            / self.T
        )
        logl_bias_to_mean = -((rates.mean() - nu_mean) ** 2) / (2 * nu_sigma**2)

        return logl_bias_to_mean * bias_to_mean
        # logl[i, a] += bias_to_mean * logl_bias_to_mean

    def penalty_expected_max_bias(
        self, params, p_N_AP, N_AP, max_spike_count, bias_to_expected_max, idx
    ):
        """
        This is not a normalized distribution!!
        """

        ### now, get extreme value distribution
        p_N_AP_cum = np.pad(
            np.cumsum(p_N_AP) ** self.data["n_neurons"][idx],
            (1, 0),
            mode="constant",
            constant_values=0,
        )
        p_extreme = np.diff(p_N_AP_cum)
        p_extreme *= 1 - 1.0 / (params["distr"][0]["nu_max"] - N_AP / self.T + 1) ** 2

        p_extreme /= p_extreme.sum()  # normalizing

        logl_extreme_empirical = np.log(p_extreme[max_spike_count])
        # print(logl_extreme_empirical)
        # print(logl_extreme_empirical)
        return bias_to_expected_max * logl_extreme_empirical

        # ## get all spike counts higher than allowed by model
        # idx_too_high = np.where(self.binom[animal]['N_AP'] > max_spike_count)[0]
        # logp_penalty = 0
        # for idx in idx_too_high:
        #     N_AP_too_high = self.binom[animal]['N_AP'][idx]
        #     k_AP_too_high = self.binom[animal]['k_AP'][idx]
        #     logp_penalty += -10**2 * k_AP_too_high*(N_AP_too_high - max_spike_count)**2
        # return logp_penalty

    def get_params_from_p(self, p_in, idx_chain=None, idx=None, biological=False):

        params = {"distr": []}

        ## this does not include structure from distribution, etc
        params_tmp = super().get_params_from_p(p_in, idx_chain=idx_chain, idx=idx)

        if biological:
            # QUESTION:
            #     * how to do this quickly (solving selfcon, ...)?
            #     * is approximation good enough?

            params["distr"][0]["nu_max"] = get_nu_max(
                params_tmp["nu_bar"],
                params_tmp["tau_A"],
                params_tmp["tau_N"],
                params_tmp["r_N"],
            )
            params["distr"][0]["gamma"] = get_gamma(
                params_tmp["nu_bar"],
                params_tmp["alpha_0"],
                params_tmp["tau_A"],
                params_tmp["tau_N"],
                params_tmp["r_N"],
            )
            params["distr"][0]["delta"] = get_delta(
                params_tmp["nu_bar"],
                params_tmp["alpha_0"],
                params_tmp["tau_A"],
                params_tmp["tau_N"],
                params_tmp["r_N"],
            )
            params = build_distr_structure_from_params(params, self.paramIn)

        else:

            params = build_distr_structure_from_params(params_tmp, self.paramIn)

            # print("params_tmp: ", params_tmp)
        # print("\n params:", params)

        return params


def build_distr_structure_from_params(params_tmp, paramNames):

    params = {"distr": []}
    for var in paramNames:
        var_split = var.split("_")
        var_root = "_".join(var_split[:-1])
        idx = int(var_split[-1][-1]) if var_split[-1].startswith("pop") else np.nan
        if len(params["distr"]) < idx + 1:
            params["distr"].append({})

        if np.isnan(idx):
            params[var] = params_tmp[var]
        else:
            params["distr"][int(idx)][var_root] = params_tmp[var]

    return params


def expected_maximum_value(p_N_AP, n, max_spike_count=None):
    """
    p_N_AP
        - probability of observing N_AP spikes in a neuron (N_AP being the index of the array)
    n
        - the number of neurons in the population

    output:
        - expected maximum value probability distribution
    """

    ### now, get extreme value distribution
    p_N_AP_cum = np.pad(
        np.cumsum(p_N_AP) ** n,
        (1, 0),
        mode="constant",
        constant_values=0,
    )
    p_extreme = np.diff(p_N_AP_cum)
    # p_extreme *= 1 - 1.0 / (params["distr"][0]["nu_max"] - N_AP / self.T + 1) ** 2
    p_extreme /= p_extreme.sum()  # normalizing

    if max_spike_count is None:
        return p_extreme
    else:
        return p_extreme[max_spike_count]

    # logl_extreme_empirical = np.log(p_extreme[max_spike_count])
    # print(logl_extreme_empirical)
    # return logl_extreme_empirical
    # return p_extreme


def weight(x, mode="sigmoid", **kwargs):
    # return a + (1-a)/np.log(np.exp(1)+1/(1-x))
    if mode == "sigmoid":
        return kwargs["offset"] + (1 - kwargs["offset"]) / (
            1 + np.exp(kwargs["slope"] * (x - kwargs["threshold"]))
        )
    else:
        return None


# def weight(x, a=0.5, b=10):
#     return a + (1 - a) * (1 - x) ** b


def get_p_N_AP(nu, p, T, correct_N=5, correct_threshold=10 ** (-4), _print=False):

    p_N_AP = p_nu(nu, p) / T
    p_N_AP[~np.isfinite(p_N_AP)] = np.nan

    if correct_N == 0:
        return p_N_AP

    ## correcting terms until normalization becomes decent
    p_full = p_N_AP[:correct_N].copy()
    for N in range(correct_N):
        ### correcting one at a time
        p_N_AP[N] = adaptive_integration(f, 0, 100.0 / T, args=(p, N, T), eps_pow=-2)

        idx_start = max(0, N - 3)
        dp = np.abs(p_full[idx_start : N + 1] - p_N_AP[idx_start : N + 1])

        ## ensure that p_N_AP converges and is not just oscillating
        if np.all(dp < correct_threshold):  # all recent corrections have been small
            # print(f"approximation good enough @{N=}, {dp[-1]}")
            break

    return p_N_AP


def get_default_priors():

    halfnorm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(x)
    norm_ppf = lambda x, mean, sigma: mean + sigma * np.sqrt(2) * erfinv(2 * x - 1)

    prior = {}
    # prior["gamma_0"] = BM.prior_structure(
    #     norm_ppf,
    #     mean=BM.prior_structure(halfnorm_ppf, loc=1.0, scale=1.0),
    #     sigma=BM.prior_structure(halfnorm_ppf, loc=0.0, scale=1.0),
    # )

    prior["gamma_0"] = prior_structure(norm_ppf, mean=2.0, sigma=0.5)
    prior["delta_0"] = prior_structure(norm_ppf, mean=6.0, sigma=2.0)
    prior["nu_max_0"] = prior_structure(
        None,
        value=25.0,
    )
    return prior


def run_sampling(
    event_counts,
    T,
    priors,
    mode="ultranest",
    key=None,
    biological=False,
    correct_N=5,
    bias_to_mean=0,
    bias_to_expected_max=0,
    n_live=100,
    nP=1,
    logLevel=logging.ERROR,
):

    # withZeros = True
    # two_pop = mP.two_pop
    two_pop = False

    BM = BayesModel(logLevel=logLevel)
    BM.prepare_data(event_counts, T)
    BM.two_pop = two_pop

    if priors is None:
        priors = get_default_priors()
    BM.set_priors(priors)

    pprint.pprint(BM.priors)
    vectorized = mode == "ultranest"
    my_prior_transform = BM.set_prior_transform(vectorized=vectorized)
    my_likelihood = BM.set_logl(
        vectorized=vectorized,
        correct_N=correct_N,
        bias_to_expected_max=bias_to_expected_max,
        bias_to_mean=bias_to_mean,
        biological=biological,
    )

    if mode=='dynesty':

        from dynesty import NestedSampler, pool as dypool

        print('running nested sampling')
        # print(np.where(hbm.pTC['wrap'])[0])
        if nP>1:
            with dypool.Pool(nP,my_likelihood,my_prior_transform) as pool:
                # sampler = DynamicNestedSampler(pool.loglike,pool.prior_transform,BM.nParams,
                #         pool=pool,
                #         # periodic=np.where(BM.wrap)[0],
                #         sample='slice'
                #     )
                # print('idx of p: ',BM.priors['p']['idx'])
                sampler = NestedSampler(
                    pool.loglike,
                    pool.prior_transform,
                    BM.n_params,
                    pool=pool,
                    nlive=n_live,
                    bound="single",
                    reflective=[BM.priors["p"]["idx"]] if two_pop else False,
                    # periodic=np.where(BM.wrap)[0],
                    sample="rslice",
                )
                sampler.run_nested(dlogz=1.)
        else:

            # import ultranest.plot as ultraplot
            # from ultranest.popstepsampler import PopulationSliceSampler, generate_region_oriented_direction

            sampler = NestedSampler(
                my_likelihood,
                my_prior_transform,
                BM.nParams,
                nlive=n_live,
                bound="single",
                reflective=[BM.priors["p"]["idx"]] if two_pop else False,
                # periodic=np.where(BM.wrap)[0],
                sample="rslice",
            )
            sampler.run_nested(dlogz=1.)
        sampling_result = sampler.results

        return BM, sampling_result, sampler
    else:
        import ultranest
        import ultranest.stepsampler

        from ultranest.popstepsampler import (
            PopulationSliceSampler,
            generate_region_oriented_direction,
        )
        from ultranest.mlfriends import RobustEllipsoidRegion

        NS_parameters = {
            "min_num_live_points": n_live,
            "max_num_improvement_loops": 3,
            "max_iters": 50000,
            "cluster_num_live_points": 20,
        }

        sampler = ultranest.ReactiveNestedSampler(
            BM.paramNames, 
            my_likelihood, my_prior_transform,
            wrapped_params=BM.wrap,
            vectorized=True,num_bootstraps=20,
            ndraw_min=512
        )

        logger = logging.getLogger("ultranest")
        logger.setLevel(logging.ERROR)

        show_status = True
        n_steps = 10
        sampler.stepsampler = PopulationSliceSampler(
            popsize=2**4,
            nsteps=n_steps,
            generate_direction=generate_region_oriented_direction,
        )

        sampling_result = sampler.run(
            **NS_parameters,
            region_class=RobustEllipsoidRegion,
            update_interval_volume_fraction=0.01,
            show_status=show_status,
            viz_callback="auto",
        )

        # num_samples = BM.nParams*100
        # num_samples = np.maximum(400,BM.nParams*100)
        # num_samples = nLive

        # sampling_result = sampler.run(
        #     min_num_live_points=num_samples,
        #     max_iters=20000,cluster_num_live_points=20,max_num_improvement_loops=3,
        #     show_status=True,viz_callback='auto')

        # plt.figure()
        # cornerplot(sampling_result)
        # plt.show(block=False)

        # plt.figure()
        # sampler.plot_trace()
        # plt.show(block=False)

        return BM, sampling_result, sampler


def compare_results(BM, sampler, mP, mode="ultranest", biological=False):

    print('data in:')
    for key in mP.params.keys():
        print(f'{key} = {mP.params[key]}')

    paraKeys = []

    if biological:
        paraKeys.extend(["nu_bar", "alpha_0", "tau_A", "tau_N", "r_N"])
        truth_values = None
    else:
        if mP.two_pop:
            paraKeys.extend("p")

        for i in range(2 if mP.two_pop else 1):
            paraKeys.extend([f"gamma_{i}", f"delta_{i}", f"nu_max_{i}"])
        # paraKeys.extend(["gamma", "delta_1", "nu_max_1"])
        # if mP.two_pop:
        #     paraKeys.extend(["gamma_2", "delta_2"])

        truth_values = []
        for distributions in mP.params["distr"]:
            # for key in paraKeys:
            for key in distributions:
                # if key in distributions.keys():
                truth_values.append(distributions[key])
            # break
            # truth_values.append(mP.params[key])

    mean = get_mean_from_sampler(sampler, paraKeys, mode=mode)

    if mode=='dynesty':

        from dynesty import utils as dyfunc, plotting as dyplot
        #

        # samples,weights = sampler.results.samples, sampler.results.importance_weights()
        # mean,cov = dyfunc.mean_and_cov(samples,weights)

        dyplot.traceplot(
            sampler.results,
            # truths=truth_values,
            truth_color="black",
            show_titles=True,
            trace_cmap="viridis",
        )
        plt.show(block=False)

        dyplot.cornerplot(
            sampler.results,
            color="dodgerblue",
            # truths=truth_values,
            show_titles=True,
        )
        plt.show(block=False)

        print("\nresults", mean)
    else:
        # mean = sampler.results['posterior']['mean']

        # from ultranest.plot import traceplot,cornerplot
        sampler.plot_trace()
        sampler.plot_corner()

        print('\nresults')
        for i,key in enumerate(sampler.results['paramnames']):
            print(
                f"{key} = {sampler.results['posterior']['mean'][i]} \pm {sampler.results['posterior']['stdev'][i]}"
            )
        plt.show(block=False)

    mP.plot_rates(param_in=mP.params)

    if biological:
        dictfilt = lambda x, y: dict([(i, x[i]) for i in x if i in set(y)])

        distribution_mean = {}
        distribution_mean["nu_max"] = get_nu_max(
            mean["nu_bar"],
            mean["tau_A"],
            mean["tau_N"],
            mean["r_N"],
        )
        distribution_mean["gamma"] = get_gamma(
            mean["nu_bar"],
            mean["alpha_0"],
            mean["tau_A"],
            mean["tau_N"],
            mean["r_N"],
        )
        distribution_mean["delta"] = get_delta(
            mean["nu_bar"],
            mean["alpha_0"],
            mean["tau_A"],
            mean["tau_N"],
            mean["r_N"],
        )
        nu_bar_in = get_nu_bar(
            gamma=mP.params["distr"][0]["gamma"],
            delta=mP.params["distr"][0]["delta"],
            nu_max=mP.params["distr"][0]["nu_max"],
        )
        nu_bar_out = get_nu_bar(
            gamma=distribution_mean["gamma"],
            delta=distribution_mean["delta"],
            nu_max=distribution_mean["nu_max"],
        )

        tau_I_in = get_tau_I(nu_max=mP.params["distr"][0]["nu_max"])
        tau_I_out = get_tau_I(nu_max=distribution_mean["nu_max"])

        alpha_0_in = get_alpha_0(
            gamma=mP.params["distr"][0]["gamma"],
            delta=mP.params["distr"][0]["delta"],
            nu_max=mP.params["distr"][0]["nu_max"],
        )
        alpha_0_out = get_alpha_0(
            gamma=distribution_mean["gamma"],
            delta=distribution_mean["delta"],
            nu_max=distribution_mean["nu_max"],
        )

        print(f"{nu_bar_in=}, {nu_bar_out=}")
        print(f"{tau_I_in=}, {tau_I_out=}")
        print(f"{alpha_0_in=}, {alpha_0_out=}")

    results_inferred = {"distr": [{}] * (2 if mP.two_pop else 1)}
    for key in mean.keys():
        print(f"{key} = {mean[key]}")

        if key == "p":
            results_inferred[key] = mean[key]
            continue

        # if key.startswith("nu"):
        #     var = key
        # else:
        key_split = key.split("_")
        idx = int(key_split[-1]) if len(key_split) > 1 else np.nan
        var = ("_").join(key_split[:-1]) if np.isfinite(idx) else key

        print(var, idx)
        results_inferred["distr"][idx][var] = mean[key]

    return results_inferred


def get_mean_from_sampler(results, paramNames, mode="ultranest", output="dict"):

    mean = {} if output == "dict" else []
    for i, key in enumerate(paramNames):
        if mode == "dynesty":
            samp = results.samples[:, i]
            weights = results.importance_weights()
        else:
            samp = results["weighted_samples"]["points"][:, i]
            weights = results["weighted_samples"]["weights"]

        if output == "dict":
            mean[key] = (samp * weights).sum()
        elif output == "list":
            mean.append((samp * weights).sum())
        # print(f"{key} mean: {mean[key]:.3f}")
    return mean


def log_binomial(n, k):
    return gammaln(n + 1) - gammaln(n - k + 1) - gammaln(k + 1)
