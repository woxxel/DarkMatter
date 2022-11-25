import numpy as np

def p_nu(NU,gamma,delta,nu_max):
    scaled_NU = np.log(NU / nu_max)
    return gamma / ( nu_max * np.sqrt( -np.pi * scaled_NU ) ) * \
        np.exp( - delta**2/2.) * ( NU / nu_max )**(gamma**2 - 1) * \
        np.cosh( gamma * delta * np.sqrt( -2 * scaled_NU ) )

def logp_nu(NU,gamma,delta,nu_max):
    scaled_NU = np.log(NU / nu_max)
    return - np.log( nu_max / gamma * np.sqrt( -np.pi * scaled_NU ) ) - delta**2 / 2 + \
        ( gamma**2 - 1 ) * scaled_NU + \
        np.log( np.cosh( gamma * delta * np.sqrt( -2 * scaled_NU ) ) )

def get_nu_bar(gamma,delta,nu_max):
    return nu_max * gamma / np.sqrt(gamma**2 + 1) * np.exp(-delta**2 / (2*(1+gamma**2)))

def get_nu_peak(gamma,delta,nu_max):
    return nu_max * np.exp( - (gamma**2 * delta**2 - 2*(gamma**2 - 1) + gamma * delta *np.sqrt(gamma**2 * delta**2 - 4*(gamma**2 - 1))) / (4 * (gamma**2 - 1)**2))

def get_chi(gamma,delta,nu_max):
    return -np.log10(get_nu_peak(gamma,delta,nu_max)/get_nu_bar(gamma,delta,nu_max))
