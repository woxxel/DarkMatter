import numpy as np

"""
    Transformation from biological parameters to distribution parameters
"""
TAU_A_DEFAULT = 0.01
TAU_N_DEFAULT = 0.2
R_N_DEFAULT = 0
TAU_M_DEFAULT = 0.01
J_0_DEFAULT = -1.0

def q_analytical(nu_bar, alpha_0, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):
    """
        this is an analytic approximation to the transcendent selfcon. equation
    """

    J_0 = J_0 * tau_M
    q_0 = alpha_0**2 / J_0**2
    eps = 1.0 / f_tau(tau_A, tau_N, r_N) * (nu_bar + q_0 / nu_bar)

    return (
        1
        / np.sqrt(1 + 2 * eps)
        * (1 + eps) ** ((1 + eps) / (1 + 2 * eps))
        * nu_bar ** ((2 + 2 * eps) / (1 + 2 * eps))
        * get_nu_max(nu_bar, tau_A, tau_N, r_N) ** (2 * eps / (1 + 2 * eps))
    )

def q_exact(nu_bar,alpha_0,tau_A=TAU_A_DEFAULT,tau_N=TAU_N_DEFAULT,r_N=R_N_DEFAULT,tau_M=TAU_M_DEFAULT,J_0=J_0_DEFAULT):

    J_0 = J_0 * tau_M
    q_0 = alpha_0**2 / J_0**2


    # eps = 1.0 / f_tau(tau_A, tau_N, r_N) * (nu_bar + q_0 / nu_bar)

    # return (
    #     1
    #     / np.sqrt(1 + 2 * eps)
    #     * (1 + eps) ** ((1 + eps) / (1 + 2 * eps))
    #     * nu_bar ** ((2 + 2 * eps) / (1 + 2 * eps))
    #     * get_nu_max(nu_bar, tau_A, tau_N, r_N) ** (2 * eps / (1 + 2 * eps))
    # )


def I_squared_nu(self, nu, q, p):
    nu = np.tile(nu,(2,1)) if np.isscalar(nu) else nu

    return - ( self.alpha(q,p)**2 + self.sigma_V(nu,p)**2 ) * np.log( (nu[p,...]/self.rate_max(nu,p))**2 * (1 + (self.alpha(q,p) / self.sigma_V(nu,p))**2) )

def I_squared_q(self, nu, q, p):
    q = self.broadcast_q(q)

    return -( self.alpha(q,p)**2 + 1./2 * self.sigma_V(nu, p)**2 ) * np.log( ( q[p,...]/self.rate_max(nu, p)**2 )**2 * (1 + 2*(self.alpha(q,p) / self.sigma_V(nu,p))**2) )


def sigma_V(nu_bar, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):

    J_0 = J_0 * tau_M
    sigma_V_A_sq = J_0**2 * nu_bar * f_tau_partial(tau_A, tau_N, 1.0 - r_N)
    sigma_V_N_sq = J_0**2 * nu_bar * f_tau_partial(tau_N, tau_A, r_N)

    # print(f"{sigma_V_A_sq=}, {sigma_V_N_sq=}")

    return np.sqrt(sigma_V_A_sq + sigma_V_N_sq)


def sigma_V_dot(nu_bar, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):

    J_0 = J_0 * tau_M
    sigma_V_A_dot_sq = (
        J_0**2
        * nu_bar
        / (tau_A * tau_M)
        * f_tau_partial(tau_A, tau_N, 1.0 - r_N, tau_M)
    )
    sigma_V_N_dot_sq = (
        J_0**2 * nu_bar / (tau_N * tau_M) * f_tau_partial(tau_N, tau_A, r_N, tau_M)
    )

    return np.sqrt(sigma_V_A_dot_sq + sigma_V_N_dot_sq)


def alpha(nu_bar, alpha_0, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):

    J_0 = J_0 * tau_M
    q = q_analytical(nu_bar, alpha_0, tau_A, tau_N, r_N)

    return np.sqrt(J_0**2 * q + alpha_0**2)


def f_tau(tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT):
    return f_tau_partial(tau_A, tau_N, 1 - r_N, tau_M) + f_tau_partial(
        tau_N, tau_A, r_N, tau_M
    )

def f_tau_partial(tau_I, tau_I_second, r, tau_M=TAU_M_DEFAULT):
    # print(f"{tau_I=}, {tau_I_second=}, {r=}")
    return (
        1.0
        / (tau_I + tau_M)
        * (r**2 / 2.0 + (r * (1.0 - r) * tau_I) / (tau_I + tau_I_second))
    )


"""
    obtain the distribution parameters
"""
def get_nu_max(nu_bar, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):
    """
    Getting the maximum firing rate value
    """
    return (
        1.0
        / (2 * np.pi)
        * sigma_V_dot(nu_bar, tau_A, tau_N, r_N, tau_M, J_0)
        / sigma_V(nu_bar, tau_A, tau_N, r_N, tau_M, J_0)
    )


def get_gamma(nu_bar, alpha_0, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):

    return sigma_V(nu_bar, tau_A, tau_N, r_N, tau_M, J_0) / alpha(
        nu_bar, alpha_0, tau_A, tau_N, r_N, tau_M, J_0
    )


def get_delta(nu_bar, alpha_0, tau_A=TAU_A_DEFAULT, tau_N=TAU_N_DEFAULT, r_N=R_N_DEFAULT, tau_M=TAU_M_DEFAULT, J_0=J_0_DEFAULT):

    gamma = get_gamma(nu_bar, alpha_0, tau_A, tau_N, r_N, tau_M, J_0)
    nu_max = get_nu_max(nu_bar, tau_A, tau_N, r_N, tau_M, J_0)

    delta_sq = -(1 + gamma**2) * np.log((nu_bar / nu_max) ** 2 * (1 / gamma**2 + 1))

    return np.sqrt(delta_sq)
