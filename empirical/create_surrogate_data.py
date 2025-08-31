
import numpy as np
from typing import List, Dict, Tuple


from scipy.integrate import quad

from matplotlib import pyplot as plt
# from DM_theory.transform_meta_to_bio import get_nu_bar,get_tau_I,get_alpha_0

# from inference.utils.utils import p_nu
# from DM_theory.network import Network
# from inference.transform_meta_to_bio import get_nu_bar,get_tau_I,get_alpha_0
from inference.utils.structures import DistributionModelParams as distr_params
from inference.BayesModel import rho_nu


class SurrogateData:

    def __init__(self,N=100,T=1200,n_animals=1):

        self.N = N
        self.T = T
        self.n_animals = n_animals


    def artificial_data(
            self,
            parameter: Dict,
            plot: bool = False,
        ):

        self.parameter = parameter

        default = [(j,) for j in range(len(parameter["distr"]))]
        join = parameter.get("join",default)
        if not join:
            join = default

        # join=parameter.get("join",[])
        # join = [] if not join else join
        rho = rho_nu(np.linspace(0,1,10),parameter["distr"],join)

        self.n_pop = 1 if len(rho.shape)==1 else rho.shape[0]

        # error_status = calculate_biological_parameters_from_meta(params=parameter["distr"])
        # if error_status:
        #     return

        self.event_counts = np.full((self.n_animals, self.n_pop, int(self.N)), np.nan)
        for a in range(self.n_animals):
            for n in range(self.n_pop):
                # print(join_)
                samples, samples_T = self.draw_samples(
                    f_target=rho_nu,
                    args=([parameter["distr"][j] for j in join[n]],True),
                    high=np.max([parameter["distr"][j].nu_max for j in join[n]]),
                    tolerance=0.001,
                    plot=plot)

                self.event_counts[a,n,:] = samples_T

        # rates = event_counts/self.T

        # if plot:
        #     self.plot_rates(key=f"{parameter['distr'][0]['gamma']=}")
        return self.event_counts


    def draw_samples(self,f_target: callable, args: Tuple, high: float, tolerance: float=0.001,plot=False,save=False):

        ## finding the lower bound for the bounding function by requiring the CDF(g(nu),0,b) < 0.001
        ## necessary, as function has a pole at nu->0 and requires starting at non-zero point
        i = -1
        while True:
            low = 10**i

            res,err = quad(f_target,0,low,args=args,points=np.logspace(-3,0,4))
            if res < tolerance:
                break

            i -= 1
        # print(f'lower bound: {low=}')

        # define M such that M*g(nu) > f(nu) for all nu
        nu = np.logspace(np.log10(low),np.log10(high),10**3+1)

        M = np.ceil(np.nanmax(f_target(nu,*args) / normalized_envelope_density(nu,low,high)))

        samples = rejection_sampling(f_target,args,self.N, M, low, high)
        # print("samples:",samples.shape,samples.max(),samples)
        samples_T = np.random.poisson(samples*self.T,samples.shape)
        # print("samples_T:",samples_T.shape,samples_T.max(),samples_T)

        if plot:
            nbins = max(101,self.N//20)
            fig,ax = plt.subplots(1,2,figsize=(6,3))
            ax[0].hist(samples,bins=np.logspace(np.log10(low),np.log10(high),nbins),density = True,label='samples',color="tab:blue",alpha=0.5)
            ax[0].hist(samples_T/self.T,bins=np.logspace(np.log10(low),np.log10(high),nbins),density = True,label='samples(T)',color="tab:red",alpha=0.5)
            ax[0].plot(nu, f_target(nu,*args),'k-',label='target f')
            ax[0].plot(nu, M*normalized_envelope_density(nu,low,high),'r-',label='envelope g')
            ax[0].axvline(low,color='r',linewidth=0.5,linestyle='--',label='lower bound')
            ax[0].axvline(1/self.T,color='k',linewidth=0.5,linestyle='--',label='min count')
            ax[0].set_xscale('log')
            ax[0].set_yscale('log')
            ax[0].legend(bbox_to_anchor=(1.05, 1.05), loc='upper right',ncol=2)

            ax[1].hist(samples,bins=np.linspace(0,high,nbins),density = True,label='samples',color="tab:blue",alpha=0.5)
            ax[1].hist(samples_T/self.T,bins=np.linspace(0,high,nbins),density = True,label='samples(T)',color="tab:red",alpha=0.5)
            ax[1].plot(nu, f_target(nu,*args),'k-')
            ax[1].plot(nu, M*normalized_envelope_density(nu,low,high),'r-')

            rho = f_target(nu,*args)
            y_max = np.max(rho[np.isfinite(rho)]) * 1.3
            # y_max = np.max(f_target(nu,*args))*1.1
            plt.setp(ax[1],xlim=[0,high/2],ylim=[0,y_max])

            for axx in ax:
                plt.setp(axx,xlabel='rate [Hz]',ylabel='density')
                axx.spines[['top','right']].set_visible(False)

            plt.tight_layout()
            if save:
                plt.savefig(f'./figures/rejection_sampling_example.png')
            plt.show(block=False)

        return samples, samples_T


    def plot_rates(self,param_in=None):

        fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))

        # gamma = gamma if gamma else self.params['gamma']
        # delta = delta if delta else self.params['delta']
        # nu_max = nu_max if nu_max else self.params['nu_max']
        # parameters = [self.params]
        # if param_in:None
        # parameters.append(param_in)

        # n_pop = self.event_counts.shape[1]
        # param=param_in
        param = param_in if param_in else self.parameter

        # rates = self.rates[key] if key else self.rates
        rates = self.event_counts/self.T
        # rates = self.rates

        xlim = 25.
        nbins = 101
        n_steps = 10001

        print(f"{self.n_pop=}")
        colors = ['k','tab:red']
        if param:
            # bins = np.linspace(0, param["distr"][0]["nu_max"], 101)

            NU_log = np.logspace(-20, np.log10(xlim), n_steps)
            # print(NU_log)
            p_NU = rho_nu(NU_log, param["distr"],param["join"])
            p_NU_cum = np.nancumsum(rho_nu(NU_log, param["distr"],param["join"])[...,:-1] * np.diff(NU_log),axis=-1)

            ax[0].plot(
                NU_log, p_NU.T, color="k", linestyle="--", label="original distribution"
            )

            for n in range(self.n_pop):
                ax[1].plot(
                    NU_log[:-1],
                    p_NU_cum[n,:],
                    label=f"original (pop {n})",
                    color=colors[n],
                    linestyle="--",
                    linewidth=2.
                )
            ax[1].legend()
        # else:
            # bins = np.linspace(0, 20, 51)
        

        for a in range(self.n_animals):
            ax[0].hist(rates[a,...].T,bins=nbins,range=(0,xlim),density=True,alpha=0.8)
            ax[1].hist(rates[a,...].T,bins=nbins,range=(0,xlim),density=True,cumulative=True,histtype='step',alpha=0.6,linewidth=0.5,color=colors)
        # ## plot histogram of empirical or artificial rates
        # ax[0].hist(rates[0,...].T,bins=nbins,range=(0,xlim),density=True)
        # # ax[0].set_xscale('log')

        # ## plot underlying original distribution
        # ax[1].hist(rates[0,...].T,bins=nbins,range=(0,xlim),density=True,cumulative=True,histtype='step')
        # # print(parameters)
        # # NU = np.linspace(0,xlim,10**8+1)
        # # for i,param in enumerate(parameters):
        # # print(p_NU)

        # # bins = 10**np.linspace(-4,2,101)
        # # p_NU[-1] = 0
        # # p_NU_cum = np.nancumsum(p_NU)
        # # p_NU_cum /= np.nanmax(p_NU_cum)
        # # print(p_NU,p_NU_cum)

        plt.setp(ax[0],xlim=[10**(-4),xlim],ylim=[0,2])
        plt.setp(ax[1],xlim=[10**(-4),xlim],ylim=[0,1.1])

        ax[0].spines[["top", "right"]].set_visible(False)
        ax[1].spines[["top", "right"]].set_visible(False)
        plt.show(block=False)


def normalized_envelope_density(x, low, high):
    return 1 / (x * np.log(high / low))

def sample_from_normalized_envelope(low, high):
    u = np.random.uniform(0, 1)
    return low * (high / low)**u  # Sample from g(x) normalized over [b, a]

def rejection_sampling(f_target,args,N, M, low, high):
    samples = []
    while len(samples) < N:
        nu = sample_from_normalized_envelope(low, high)
        v = np.random.uniform(0, 1)
        if v <= f_target(nu,*args) / (M * normalized_envelope_density(nu, low, high)):
            samples.append(nu)
    return np.array(samples)


# def calculate_biological_parameters_from_meta(params):
#     """
#         returns True, if any parameter is NaN
#     """

#     if not isinstance(params,list):
#         params = [params]
    
#     bio_params = {
#         "tau_M": 0.01,
#         "J0": [-1.],
#         "rateWnt": [],
#         "tau_I": [],
#         "alpha_0": []
#     }

#     if len(params)>1:
#         print("Theoretical inference from multiple distributions not possible")
#         # return False


#     for m,param in enumerate(params):
#         bio_params["rateWnt"].append(
#             get_nu_bar(param)
#         )
#         if len(params)==1:
#             bio_params["tau_I"].append(
#                 get_tau_I(
#                     nu_max=param.nu_max, tau_m=bio_params["tau_M"]
#                 )
#             )
#             bio_params["alpha_0"].append(
#                 get_alpha_0(
#                     param,
#                     tau_m=bio_params["tau_M"],
#                     J_0=bio_params["J0"][0],
#                 )
#             )
#     # options['mode_selfcon'] = 0
#     # print(options)
#     # print('\n')
#     string_input_params = "input parameters: "
#     for distribution in params:
#         for key,val in distribution.__dict__.items():
#             string_input_params += f"{key}={val}, "
#     print(string_input_params)

#     string_inferred_params = "inferred parameters: "
#     for key in ["rateWnt", "tau_I", "alpha_0"]:
#         string_inferred_params += f"{key}={bio_params[key]}, "
#     print(string_inferred_params)

#     if np.isnan(bio_params['rateWnt']).any() or np.isnan(bio_params['tau_I']).any() or np.isnan(bio_params['alpha_0']).any():
#         print('inferred parameters contain NaNs - breaking!')
#         # self.rates = False
#         return True
#     return False
