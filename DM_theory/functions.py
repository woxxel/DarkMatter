import numpy as np

## defining functions
def p_nu(NU,gamma,delta,nu_max):

    return gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) * \
        np.exp( - delta**2/2.) * ( NU / nu_max )**(gamma**2 - 1) * \
        np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max) ) )

def get_nu_peak(gamma,delta,nu_max):
    return nu_max * np.exp( - (gamma**2 * delta**2 - 2*(gamma**2 - 1) + gamma * delta *np.sqrt(gamma**2 * delta**2 - 4*(gamma**2 - 1))) / (4 * (gamma**2 - 1)**2))

def get_chi(gamma,delta,nu_max):
    return -np.log10(get_nu_peak(gamma,delta,nu_max)/get_nu_bar(gamma,delta,nu_max))

def get_tau_I(nu_max,tau_m=0.01):
    return 1./tau_m * (1./ (2 * np.pi * nu_max))**2

def get_nu_bar(gamma,delta,nu_max):
    return nu_max * gamma / np.sqrt(gamma**2 + 1) * np.exp(-delta**2 / (2*(1+gamma**2)))

def get_q(gamma,delta,nu_max):
    return nu_max**2 * gamma / np.sqrt(gamma**2 + 2) * np.exp(-delta**2 / (2+gamma**2))

def get_alpha_0(gamma,delta,nu_max,tau_m=0.01,J_0=-1):
    
    J_0 *= tau_m
    
    nu_mean = get_nu_bar(gamma,delta,nu_max)
    tau_I = get_tau_I(nu_max,tau_m)

    q = get_q(gamma,delta,nu_max)
    # print(f'{nu_mean=}, {tau_I=}, {q=}')
    # print(nu_mean,tau_I,q)
    # print(tau_I * nu_mean/ (2 * gamma**2 * (tau_I + tau_m)))
    
    return np.sqrt(J_0**2 * ( nu_mean/ (2 * gamma**2 * (tau_I + tau_m)) - q))


def border1(gamma,delta,nu_max):

    nu_mean = get_nu_bar(gamma,delta,nu_max)
    return (nu_mean/nu_max)**2 * (1/gamma**2 + 1)

def border2(gamma,delta,nu_max):

    q = get_q(gamma,delta,nu_max)
    return (q / nu_max**2)**2 * (2/gamma**2 + 1)