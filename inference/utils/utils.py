import numpy as np

from scipy.integrate import quad
from scipy.special import gammaln

from .structures import DistributionModelParams as distr_params


# def p_nu_single(NU,gamma,delta,nu_max,log=False):

#     if not isinstance(NU,np.ndarray):
#         NU = np.array(NU)

#     if log:
#         # return - np.log( nu_max / gamma * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) - delta**2 / 2 + \
#         ( gamma**2 - 1 ) * np.log( NU / nu_max ) + \
#         np.log( np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max ) ) ) )
#     else:
#         p = np.zeros_like(NU)
#         NU_mask = NU > 0
#         NU_scaled = NU[NU_mask] / nu_max
#         # NU_scaled = NU / nu_max
#         p[NU_mask] = gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU_scaled ) ) ) * \
#             np.exp( - delta**2/2.) * ( NU_scaled )**(gamma**2 - 1) * \
#             np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU_scaled ) ) )
#         p[~NU_mask] = 0 if gamma > 1 else np.inf
        
#         return p

# def p_nu(NU,p,log=False):
#     '''
#         calculates the probability to observe a firing rate NU, given the underlying distribution
#     '''

#     if 'p' in p: ## two populations
#         # return (p['weight_dark'] * p_nu_single(NU,p['gamma_dark'],p['delta_dark'],p['nu_max']) + \
#         # (1-p['weight_dark']) * p_nu_single(NU,p['gamma'],p['delta'],p['nu_max']))
#         p_nu = 0
#         prob = np.copy(p['p'])
#         for m in range(2):
#             p_nu += prob * p_nu_single(
#                 NU,
#                 p["distr"][m]["gamma"],
#                 p["distr"][m]["delta"],
#                 p["distr"][m]["nu_max"],
#                 log=log,
#             )
#             prob = 1 - prob

#         return p_nu
#     # (p['p'] * p_nu_single(NU,p['gamma'],p['delta'],p['nu_max'],log=log) + \
#     # (1-p['p']) * p_nu_single(NU,p['gamma2'],p['delta2'],p['nu_max2'],log=log))

#     else:
#         m=0
#         return p_nu_single(
#             NU,
#             p["distr"][m]["gamma"],
#             p["distr"][m]["delta"],
#             p["distr"][m]["nu_max"],
#             log=log,
#         )

def rho_nu(NU: np.ndarray, params: distr_params, join = False, log: bool=False):

    params = params if isinstance(params,list) else [params]
    n_pop = len(params)
    NU = np.atleast_1d(NU)

    rho = np.zeros((n_pop,len(NU)))

    NU_mask = NU > 0
    for i,p in enumerate(params):
        NU_scaled = NU[NU_mask] / p.nu_max

        if log:
            rho[i,NU_mask] = ( p.gamma**2 - 1 ) * np.log( NU / p.nu_max ) + \
                np.log( np.cosh( p.gamma * p.delta * np.sqrt( -2 * np.log( NU / p.nu_max ) ) ) )
            rho[i,~NU_mask] = -np.inf if p.gamma > 1 else np.inf
        else:
            rho[i,NU_mask] = p.gamma / ( p.nu_max * np.sqrt( -np.pi * np.log( NU_scaled ) ) ) * \
                np.exp( - p.delta**2/2.) * ( NU_scaled )**(p.gamma**2 - 1) * \
                np.cosh( p.gamma * p.delta * np.sqrt( -2 * np.log( NU_scaled ) ) )
            rho[i,~NU_mask] = 0 if p.gamma > 1 else np.inf

    if n_pop==1:
        return rho[0,:]
    
    if isinstance(join,bool):
        join = [tuple(range(n_pop))] if join else []

    if len(join) > 0:
        rho_out = np.zeros((len(join),len(NU)))
        for j,join_ in enumerate(join):
            kappa = np.array([params[i].kappa for i in join_])
            kappa /= kappa.sum()

            rho_out[j,:] = (kappa[:,np.newaxis] * rho[join_,:]).sum(axis=0)
        return np.squeeze(rho_out)
    else:
        return rho


def poisson_spikes(nu,N_AP,T,log=False):
    ## using the gamma-function to obtain log(N!) = gammaln(N+1)
    # nu = nu[:,np.newaxis]
    if log:
        return N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T
    else:
        return np.exp(N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T)


def spike_observation_probability(NU,args_rho,args_poisson,log=False):
    '''
        calculates the probability to observe N_AP action potentials in any neuron, 
        given the underlying firing rate distribution and poisson firing
    '''    
    
    # missing: proper handover of join, proper evaluation, ....
    if log:
        return np.exp(rho_nu(NU,*args_rho,log=log) + poisson_spikes(NU,*args_poisson,log=log))
    else:
        return rho_nu(NU,*args_rho,log=log) * poisson_spikes(NU,*args_poisson,log=log)


def adaptive_integration(f,x_lower,x_upper,args,eps_pow=-8,eps_thr=-4):
    # while True:
        # if eps_pow==eps_thr:
        #     print('tolerance too high - breaking!')
        #     return None

    # print("Starting adaptive integration...")
    # print(x_lower, x_upper, args, eps_pow)
    # try:
    res,err = quad(f,x_lower,x_upper,args=args,epsabs=10**eps_pow, epsrel=10**eps_pow,points=np.logspace(-3,0,4))
    # except:
    #     res = np.nan
        # break
        # except:
        #     print(f'error in integration with tolerance 10^{eps_pow}')
        #     eps_pow += 1
    return res#,err
