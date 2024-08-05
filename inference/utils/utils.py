import numpy as np

from scipy.integrate import quad
from scipy.special import gammaln


def p_nu_single(NU,gamma,delta,nu_max,log=False):

    if not isinstance(NU,np.ndarray):
        NU = np.array(NU)

    if log:
        # return - np.log( nu_max / gamma * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) - delta**2 / 2 + \
        ( gamma**2 - 1 ) * np.log( NU / nu_max ) + \
        np.log( np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max ) ) ) )
    else:
        p = np.zeros_like(NU)
        NU_mask = NU > 0
        NU_scaled = NU[NU_mask] / nu_max
        # NU_scaled = NU / nu_max
        p[NU_mask] = gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU_scaled ) ) ) * \
            np.exp( - delta**2/2.) * ( NU_scaled )**(gamma**2 - 1) * \
            np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU_scaled ) ) )
        p[~NU_mask] = 0 if gamma > 1 else np.inf
        
        return p

def p_nu(NU,p,two_pop=False,log=False):

    if two_pop:
        # return (p['weight_dark'] * p_nu_single(NU,p['gamma_dark'],p['delta_dark'],p['nu_max']) + \
        # (1-p['weight_dark']) * p_nu_single(NU,p['gamma'],p['delta'],p['nu_max']))
        return (p['weight_dark'] * p_nu_single(NU,p['gamma_dark'],p['delta_dark'],p['nu_max'],log=log) + \
        (1-p['weight_dark']) * p_nu_single(NU,p['gamma'],p['delta'],p['nu_max'],log=log))

    else:
        return p_nu_single(NU,p['gamma'],p['delta'],p['nu_max'],log=log)

def poisson_spikes(nu,N_AP,T,log=False):
    ## using the gamma-function to obtain log(N!) = gammaln(N+1)
    if log:
        return N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T
    else:
        return np.exp(N_AP*np.log(nu*T) - gammaln(N_AP+1) - nu*T)


def f(NU,p,N_AP,T,two_pop=False,log=False):
    '''
        calculates the probability to observe N_AP action potentials in any neuron, 
        given the underlying firing rate distribution and poisson firing
    '''
    if log:
        p = np.exp(p_nu(NU,p,two_pop=two_pop,log=log) + poisson_spikes(NU,N_AP,T,log=log))
        return p
    else:
        return p_nu(NU,p,two_pop=two_pop,log=log) * poisson_spikes(NU,N_AP,T,log=log)
    
# def f(NU,p,N_AP,T,zero=False,two_pop=False,log=False):
#     '''
#         calculates the probability to observe N_AP action potentials in any neuron, 
#         given the underlying firing rate distribution and poisson firing
#     '''
#     if log:#
#         if zero:
#             return np.exp(p_nu(NU,p,two_pop=two_pop,log=log) - NU*T)
#         else:
#             # return np.exp(p_nu(NU,p,two_pop=two_pop,log=log) + poisson_spikes(NU,N_AP[:,np.newaxis],T,log=log))
#             return np.exp(p_nu(NU,p,two_pop=two_pop,log=log) + poisson_spikes(NU,N_AP,T,log=log))
#     else:
#         if zero:
#             return p_nu(NU,p,two_pop=two_pop,log=log) * np.exp(-NU*T)
#         else:
#             # return p_nu(NU,p,two_pop=two_pop,log=log) * poisson_spikes(NU,N_AP[:,np.newaxis],T,log=log)
#             return p_nu(NU,p,two_pop=two_pop,log=log) * poisson_spikes(NU,N_AP,T,log=log)



def adaptive_integration(f,x_lower,x_upper,args,eps_pow=-8,eps_thr=-4):
    while True:
        if eps_pow==eps_thr:
            print('tolerance too high - breaking!')
            return None

        # try:
        res,err = quad(f,x_lower,x_upper,args=args,epsabs=10**eps_pow, epsrel=10**eps_pow,points=np.logspace(-3,0,4))

        break
        # except:
        #     print(f'error in integration with tolerance 10^{eps_pow}')
        #     eps_pow += 1
    return res
