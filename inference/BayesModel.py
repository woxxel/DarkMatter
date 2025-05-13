import os, logging
import numpy as np

from scipy.special import binom as sp_binom, erfinv

from .HierarchicalModelDefinition import HierarchicalModel

from collections import Counter

from matplotlib import pyplot as plt


from .utils.utils import p_nu, adaptive_integration, f

from DM_theory.functions import get_nu_bar, get_q, get_tau_I, get_alpha_0
from DM_theory.functions import get_nu_max, get_gamma, get_delta

import warnings
warnings.filterwarnings("ignore")


logging.basicConfig(level=logging.INFO)


class BayesModel(HierarchicalModel):

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

    def set_priors(self, hierarchical=[], two_pop=False, biological=False):

        '''
            builds prior distributions for the model parameters
            each parameter can either be hierarchical or not
        '''

        self.two_pop = two_pop

        ## define distributions inline, instead of using scipy ones for the sake of performance
        halfnorm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(x)
        norm_ppf = lambda x, loc, scale: loc + scale * np.sqrt(2) * erfinv(2*x - 1)

        self.priors_init = {}

        if biological:
            self.priors_init[f"nu_bar"] = {
                "hierarchical": {
                    "function": norm_ppf,
                    "params": {"loc": "mean", "scale": "sigma"},
                },
                "mean": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0.0, "scale": 1.0},
                },
                "sigma": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0, "scale": 0.1},
                },
            }
            self.priors_init[f"alpha_0"] = {
                "hierarchical": {
                    "function": norm_ppf,
                    "params": {"loc": "mean", "scale": "sigma"},
                },
                "mean": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0.0, "scale": 0.1},
                },
                "sigma": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0, "scale": 0.02},
                },
            }
            self.priors_init[f"tau_A"] = {
                "hierarchical": {
                    "function": norm_ppf,
                    "params": {"loc": "mean", "scale": "sigma"},
                },
                "mean": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0.0, "scale": 0.02},
                },
                "sigma": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0, "scale": 0.005},
                },
            }
            self.priors_init[f"tau_N"] = {
                "hierarchical": {
                    "function": norm_ppf,
                    "params": {"loc": "mean", "scale": "sigma"},
                },
                "mean": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0.0, "scale": 0.2},
                },
                "sigma": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0, "scale": 0.01},
                },
            }
            self.priors_init[f"r_N"] = {
                "hierarchical": {
                    "function": norm_ppf,
                    "params": {"loc": "mean", "scale": "sigma"},
                },
                "mean": {
                    "function": lambda x: x,
                    "params": {},
                },
                "sigma": {
                    "function": halfnorm_ppf,
                    "params": {"loc": 0, "scale": 0.1},
                },
            }

        else:
            if two_pop:
                self.priors_init["p"] = {
                    "hierarchical": {
                        "function": norm_ppf,
                        "params": {"loc": "mean", "scale": "sigma"},
                    },
                    "mean": {
                        "params": {},
                        "function": lambda x: x * 0.5,
                    },
                    "sigma": {
                        "params": {"loc": 0.0, "scale": 0.1},
                        "function": halfnorm_ppf,
                    },
                }

            for m in range(2 if two_pop else 1):
                self.priors_init[f"gamma_{m}"] = {
                    "hierarchical": {
                        "function": norm_ppf,
                        "params": {"loc": "mean", "scale": "sigma"},
                    },
                    "mean": {
                        "function": norm_ppf,
                        "params": {"loc": 2.0, "scale": 0.5},
                    },
                    "sigma": {
                        "function": halfnorm_ppf,
                        "params": {"loc": 0, "scale": 0.1},
                    },
                }
                self.priors_init[f"delta_{m}"] = {
                    "hierarchical": {
                        "function": norm_ppf,
                        "params": {"loc": "mean", "scale": "sigma"},
                    },
                    "mean": {
                        "function": norm_ppf,
                        "params": {"loc": 6, "scale": 2},
                    },
                    "sigma": {
                        "function": halfnorm_ppf,
                        "params": {"loc": 0, "scale": 0.5},
                    },
                }
                if m == 0:
                    self.priors_init[f"nu_max_{m}"] = {
                        "hierarchical": {
                            "function": norm_ppf,
                            "params": {"loc": "mean", "scale": "sigma"},
                        },
                        "mean": {
                            "function": halfnorm_ppf,
                            "params": {"loc": 1.0, "scale": 50.0},
                        },
                        "sigma": {
                            "function": halfnorm_ppf,
                            "params": {"loc": 0, "scale": 5},
                        },
                    }

        '''
            translating the above specified priors into a format that can be used by the model
        '''
        super().set_priors(self.priors_init, hierarchical=hierarchical)

    def set_logl(
        self,
        vectorized=False,
        correct_N=0,
        bias_to_expected_max=0,
        bias_to_mean=0,
        correct_threshold=0.01,
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
            max_spike_count = np.nanmax(self.spike_counts,axis=0).astype('int')

            nChain = p_in.shape[0]

            logl = np.zeros((nChain,self.nAnimals))

            for animal in range(self.nAnimals):

                rates = self.rates[np.isfinite(self.rates[:, animal]), animal]
                max_spike_count = int(np.nanmax(self.spike_counts[:, animal]))

                for i in range(nChain):

                    params = {"distr": []}

                    if biological:
                        """
                        first, calculate distribution parameters from biological parameters
                        QUESTION:
                            * how to do this quickly (solving selfcon, ...)?
                            * is approximation good enough?
                        """
                        params_tmp = {}
                        for var in self.priors_init.keys():
                            params_tmp[var] = p_in[
                                i,
                                self.priors[var]["idx"]
                                + (animal if self.priors[var]["n"] > 1 else 0),
                            ]

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

                        if (
                            np.isnan(params["distr"][0]["nu_max"])
                            or np.isnan(params["distr"][0]["gamma"])
                            or np.isnan(params["distr"][0]["delta"])
                        ):
                            logl[i, animal] = -(10**6)
                            continue
                    else:
                        for var in self.priors_init.keys():
                            var_split = var.split("_")
                            var_root = "_".join(var_split[:-1])
                            idx = int(var_split[-1])
                            if len(params["distr"]) < idx + 1:
                                params["distr"].append({})

                            if np.isfinite(idx):
                                params["distr"][int(idx)][var_root] = p_in[
                                    i,
                                    self.priors[var]["idx"]
                                    + (animal if self.priors[var]["n"] > 1 else 0),
                                ]
                            else:
                                params[var] = p_in[
                                    i,
                                    self.priors[var]["idx"]
                                    + (animal if self.priors[var]["n"] > 1 else 0),
                                ]
                    # print(f"{params_tmp=}")
                    # print(f"{p=}")

                    if self.two_pop:
                        params["distr"][1]["nu_max"] = params["distr"][0]["nu_max"]

                    # if p["distr"][0]["nu_max"] < 0:
                    #     logl[i, a] = -(10**6) * abs(p["distr"][0]["nu_max"])
                    #     continue

                    # if p['gamma'] < p['gamma2']:
                    #     logl[i,a] += -10**6 * (p['gamma2'] - p['gamma'])

                    ## calculate maximum number of spikes:
                    N_max = (params["distr"][0]["nu_max"] * self.T).astype("int")

                    N_AP = np.arange(N_max).astype("int")
                    N_AP_empirical = self.binom[animal]["N_AP"][
                        self.binom[animal]["N_AP"] < N_max
                    ]
                    k_AP_empirical = self.binom[animal]["k_AP"][
                        self.binom[animal]["N_AP"] < N_max
                    ]

                    k_AP = np.zeros(N_max, dtype="int")
                    k_AP[N_AP_empirical] = k_AP_empirical

                    ### calculate actual log-likelihood
                    p_N_AP = get_p_nu(
                        (N_AP + offset) / self.T,
                        params,
                        self.T,
                        correct_N=correct_N,
                        correct_threshold=correct_threshold,
                    )

                    binom = self.binom[animal]["binomial_coefficient"][k_AP]
                    logp = (
                        np.log(binom)
                        + k_AP * np.log(p_N_AP[N_AP])
                        + (self.data["nNeurons"][animal] - k_AP)
                        * np.log(1 - p_N_AP[N_AP])
                    )

                    logl_measures = logp.sum()
                    logl[i, animal] += logl_measures

                    ### add penalties
                    logl[i, animal] += self.penalty_nu_max(N_max, animal)

                    if bias_to_mean > 0:
                        logl[i, animal] += self.penalty_mean_bias(
                            params, bias_to_mean, animal
                        )

                    if bias_to_expected_max > 0 and (max_spike_count < N_max):
                        logl[i, animal] += self.penalty_expected_max_bias(
                            params,
                            p_N_AP,
                            N_AP,
                            max_spike_count,
                            bias_to_expected_max,
                            animal,
                        )
                    elif bias_to_expected_max > 0:
                        logl[i, animal] += bias_to_expected_max * -(10**2)

                    """
                        sum up the log likelihoods to obtain the overall estimate for this animal
                    """
                    # self.log.debug((f'logls: {logl_low=}, {logl_intermediate=}, {logl_tail=}'))

                    # logl[i,a] = logl_low + logl_intermediate + logl_tail
                    if not np.isfinite(logl[i, animal]):
                        logl[i, animal] = -(10 ** (-6))

            if vectorized:
                return logl.sum(axis=1)
            # return logl.sum(axis=1)
            else:
                self.log.debug(("logl:", logl[0, :].sum()))
                return logl[0, :].sum()

        return loglikelihood

    def penalty_nu_max(self, N_max, animal):

        ## get all spike counts higher than allowed by model
        idx_too_high = np.where(self.binom[animal]["N_AP"] > N_max)[0]
        logp_penalty = 0
        for idx in idx_too_high:
            # print(idx)
            N_AP_too_high = self.binom[animal]["N_AP"][idx]
            k_AP_too_high = self.binom[animal]["k_AP"][idx]
            # print(f'{N_AP_too_high=}, {k_AP_too_high=}')
            logp_penalty += -(10**2) * k_AP_too_high * (N_AP_too_high - N_max) ** 2
        return logp_penalty

    def penalty_mean_bias(self, params, bias_to_mean, animal):

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

        nu_sigma = np.sqrt(nu_SD - nu_mean**2) / np.sqrt(self.data["nNeurons"][animal])

        rates = self.rates[np.isfinite(self.rates[:, animal]), animal]
        logl_bias_to_mean = -((rates.mean() - nu_mean) ** 2) / (2 * nu_sigma**2)

        return logl_bias_to_mean * bias_to_mean
        # logl[i, a] += bias_to_mean * logl_bias_to_mean

    def penalty_expected_max_bias(
        self, params, p_N_AP, N_AP, max_spike_count, bias_to_expected_max, animal
    ):

        ### now, get extreme value distribution
        p_N_AP_cum = np.pad(
            np.cumsum(p_N_AP) ** self.data["nNeurons"][animal],
            (1, 0),
            mode="constant",
            constant_values=0,
        )
        p_extreme = np.diff(p_N_AP_cum)
        p_extreme *= 1 - 1.0 / (params["distr"][0]["nu_max"] - N_AP / self.T + 1) ** 2

        logl_extreme_empirical = np.log(p_extreme[max_spike_count])
        return bias_to_expected_max * logl_extreme_empirical

        # ## get all spike counts higher than allowed by model
        # idx_too_high = np.where(self.binom[animal]['N_AP'] > max_spike_count)[0]
        # logp_penalty = 0
        # for idx in idx_too_high:
        #     N_AP_too_high = self.binom[animal]['N_AP'][idx]
        #     k_AP_too_high = self.binom[animal]['k_AP'][idx]
        #     logp_penalty += -10**2 * k_AP_too_high*(N_AP_too_high - max_spike_count)**2
        # return logp_penalty


def get_p_nu(nu, p, T, correct_N=5, correct_threshold=0.01, _print=False):

    p_N_AP = p_nu(nu, p) / T
    p_N_AP[~np.isfinite(p_N_AP)] = np.nan
    # print(p_N_AP)
    ## correcting terms until normalization becomes decent
    for N in range(correct_N):
        if _print:
            print(f"{N=} , {np.nansum(p_N_AP)=}")
        if abs(np.nansum(p_N_AP) - 1.0) < correct_threshold:
            if _print:
                print("break")
            break
        p_N_AP[N] = adaptive_integration(f, 0, 100.0 / T, args=(p, N, T), eps_pow=-2)

    return p_N_AP


def run_sampling(
    mP,
    mode="ultranest",
    biological=False,
    correct_N=5,
    bias_to_mean=0,
    bias_to_expected_max=0,
    n_live=100,
    nP=1,
    logLevel=logging.ERROR,
):

    # withZeros = True
    two_pop = mP.two_pop

    BM = BayesModel(mP,mode='rates')
    # BM.set_logLevel(logging.ERROR)
    BM.set_logLevel(logLevel)
    # BM.prepare_data(mP,mode='rates',key='WT')
    BM.prepare_data(mP,mode='rates')
    BM.set_priors(hierarchical=[], biological=biological, two_pop=two_pop)
    # BM.set_priors(hierarchical=['gamma','delta','nu_max'],two_pop=two_pop)
    # BM.set_priors(hierarchical=['gamma','delta'],two_pop=two_pop)

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
                    BM.nParams,
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
        paraKeys.extend(["gamma_1", "delta_1", "nu_max_1"])
        if mP.two_pop:
            paraKeys.extend(["gamma_2", "delta_2"])

        truth_values = []
        for key in paraKeys:
            truth_values.append(mP.params[key])

    mean = {}
    for i, key in enumerate(BM.paramNames):
        if mode == "dynesty":
            samp = sampler.results.samples[:, i]
            weights = sampler.results.importance_weights()
        else:
            samp = sampler.results["weighted_samples"]["points"][:, i]
            weights = sampler.results["weighted_samples"]["weights"]

        mean[key] = (samp * weights).sum()
        print(f"{key} mean: {mean[key]}")

    if mode=='dynesty':

        from dynesty import utils as dyfunc, plotting as dyplot

        # samples,weights = sampler.results.samples, sampler.results.importance_weights()
        # mean,cov = dyfunc.mean_and_cov(samples,weights)

        dyplot.traceplot(
            sampler.results,
            truths=truth_values,
            truth_color="black",
            show_titles=True,
            trace_cmap="viridis",
        )
        plt.show(block=False)

        dyplot.cornerplot(
            sampler.results, color="dodgerblue", truths=truth_values, show_titles=True
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

    print(f'{nu_bar_in=}, {nu_bar_out=}')
    print(f'{tau_I_in=}, {tau_I_out=}')
    print(f'{alpha_0_in=}, {alpha_0_out=}')

    for key in mean.keys():
        print(f"{key} = {mean[key]}")
