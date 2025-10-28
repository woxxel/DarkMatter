import numpy as np
from .utils.structures import DistributionModelParams as distr_params

"""
    Transformation from distribution parameters to biological parameters

    inputs:
        gamma:
        delta:
        nu_max:
"""

### first (nu) and second (q) moment of firing rate
def get_nu_bar(param: distr_params):
    return param.nu_max * param.gamma / np.sqrt(param.gamma**2 + 1) * np.exp(-param.delta**2 / (2*(1+param.gamma**2)))

def get_q(param: distr_params):
    return param.nu_max**2 * param.gamma / np.sqrt(param.gamma**2 + 2) * np.exp(-param.delta**2 / (2+param.gamma**2))

### get biological parameters
def get_tau_I(param: distr_params, tau_m=0.01):
    # only holds for 1 population?!
    return 1./tau_m * (1./ (2 * np.pi * param.nu_max))**2


def get_alpha_0(param: distr_params, p=1.0, tau_m=0.01, J_0=-1.0, nP=None):

    """
        should be changed such that function manipulates parameters directly, without need to hand back
    """
    dims, kwargs = get_input_dimensions(
        nP, gamma=param.gamma, delta=param.delta, nu_max=param.nu_max, tau_m=tau_m, J_0=J_0
    )

    # print(f"{gamma=}, {delta=}, {nu_max=}, {p=}, {dims=}")
    params = distr_params(
        gamma=kwargs["gamma"], delta=kwargs["delta"], nu_max=kwargs["nu_max"]
    )
    # gamma = kwargs["gamma"]
    # delta = kwargs["delta"]
    # nu_max = kwargs["nu_max"]
    tau_m = kwargs["tau_m"]

    J_0 = kwargs["J_0"]
    J_0 *= tau_m
    # for arg in kwargs:
    #     print(arg,kwargs[arg].shape)

    nu_mean_eff = get_effective_variable(get_nu_bar, params, p, dims)
    q_eff = get_effective_variable(get_q, params, p, dims)

    tau_I = get_tau_I(params, tau_m)
    # print(f'{nu_mean_eff.shape=}, {q_eff.shape=}, {tau_I.shape=}, {gamma.shape=}')
    # print(q_eff)
    # return np.sqrt(J_0**2 * ( nu_mean_eff[...,np.newaxis]/ (2 * gamma**2 * (tau_I[...,np.newaxis] + tau_m)) - q_eff[...,np.newaxis]))
    return np.sqrt(
        J_0**2
        * (
            nu_mean_eff[..., np.newaxis] / (2 * params.gamma**2 * (tau_I + tau_m))
            - q_eff[..., np.newaxis]
        )
    )


def get_dPsi(gamma, delta, nu_max, p, tau_m=0.01, J_0=-1.0, nP=2):

    # gamma,delta,nu_max = check_arrays(gamma,delta,nu_max)
    dims, kwargs = get_input_dimensions(
        nP, gamma=gamma, delta=delta, nu_max=nu_max, tau_m=tau_m, J_0=J_0
    )
    gamma = kwargs["gamma"]
    delta = kwargs["delta"]
    nu_max = kwargs["nu_max"]
    tau_m = kwargs["tau_m"]
    J_0 = kwargs["J_0"]

    nP = dims[-1]
    # print(dims)
    assert nP == 2, "This function is only defined for 2 populations"

    J_0 *= tau_m

    # nu_mean_eff = get_effective_variable(get_nu_bar,gamma,delta,nu_max,p)
    q_eff = get_effective_variable(get_q, gamma, delta, nu_max, p, dims=dims)
    alpha_0 = get_alpha_0(gamma, delta, nu_max, p, tau_m=tau_m, J_0=J_0, nP=nP)

    # print(f'{q_eff.shape=}, {alpha_0.shape=}, {gamma.shape=}')

    Psi = delta * np.sqrt(J_0**2 * q_eff[..., np.newaxis] + alpha_0)
    # print(Psi.shape)
    return Psi[..., 1] - Psi[..., 0]


### the firing rate distribution function
def p_nu(NU,gamma,delta,nu_max):

    return gamma / ( nu_max * np.sqrt( -np.pi * np.log( NU / nu_max ) ) ) * \
        np.exp( - delta**2/2.) * ( NU / nu_max )**(gamma**2 - 1) * \
        np.cosh( gamma * delta * np.sqrt( -2 * np.log( NU / nu_max) ) )


### further parameter of interest
def get_nu_peak(gamma,delta,nu_max):
    return nu_max * np.exp( - (gamma**2 * delta**2 - 2*(gamma**2 - 1) + gamma * delta *np.sqrt(gamma**2 * delta**2 - 4*(gamma**2 - 1))) / (4 * (gamma**2 - 1)**2))

def get_chi(gamma,delta,nu_max):
    return -np.log10(get_nu_peak(gamma,delta,nu_max)/get_nu_bar(gamma,delta,nu_max))


### helper functions for calculating solutions of DM theory
def border1(gamma,delta,nu_max):

    nu_mean = get_nu_bar(gamma,delta,nu_max)
    return (nu_mean/nu_max)**2 * (1/gamma**2 + 1)

def border2(gamma,delta,nu_max):

    q = get_q(gamma,delta,nu_max)
    return (q / nu_max**2)**2 * (2/gamma**2 + 1)


# def get_alpha_0(gamma,delta,nu_max,tau_m=0.01,J_0=-1):

#     J_0 *= tau_m

#     nu_mean = get_nu_bar(gamma,delta,nu_max)
#     tau_I = get_tau_I(nu_max,tau_m)

#     q = get_q(gamma,delta,nu_max)
#     # print(f'{nu_mean=}, {tau_I=}, {q=}')
#     # print(nu_mean,tau_I,q)
#     # print(tau_I * nu_mean/ (2 * gamma**2 * (tau_I + tau_m)))

#     return np.sqrt(J_0**2 * ( nu_mean/ (2 * gamma**2 * (tau_I + tau_m)) - q))


# def get_2mod_tau_I(nu_max,tau_m=0.01):
#     return 1./tau_m * (1./ (2 * np.pi * nu_max))**2
def check_var(var, n, dims=None):

    if not hasattr(var, "__len__"):
        var = np.atleast_1d(var)

    if isinstance(var, list):
        var = np.array(var)

    # if nP is None:
    #     nP = var.shape[-1]
    nP = dims[-1]
    if len(var.shape) == 1:
        # print("test dimensions", dims, var.shape)
        assert var.shape[-1] == nP, "input is not consistent with number of populations"

        var = var[np.newaxis, :]
    else:
        if var.shape[-1] != nP:
            var = var[..., np.newaxis]

    if n == 1 and var.shape[0] > 1:
        n = var.shape[0]

    assert var.shape[0] in [1, n], "wrong dimension!"

    if var.shape[0] != n:
        var = np.tile(var, (n, 1))

    return var


def get_effective_variable(calc_fun, params, p, dims):

    n = dims[0]
    nP = dims[-1]
    if nP == 2:
        var = np.zeros((nP, n))
        for k, (g, d, n) in enumerate(
            zip(params.gamma.T, params.delta.T, params.nu_max.T)
        ):
            var[k, :] = calc_fun(g, d, n)
        return p * var[0, ...] + (1 - p) * var[1, ...]
    else:
        # return calc_fun(params.gamma[..., 0], params.delta[..., 0], params.nu_max)
        return calc_fun(params)


def check_arrays(gamma, delta, nu_max, dims=None):
    """
    sanity checks and preparation of array

    gamma and delta should be 2D arrays with
    - 1st dimension: number of data points
    - 2nd dimension: number of populations
    """

    nu_max = np.atleast_1d(nu_max)
    n = nu_max.shape[0]

    gamma = check_var(gamma, n, dims)
    delta = check_var(delta, n, dims)

    assert gamma.shape == delta.shape, "gamma and delta should have the same shape"
    # assert len(gamma.shape) == 2 and len(delta.shape) == 2, 'gamma and delta should be 2D arrays'

    return gamma, delta, nu_max


def get_input_dimensions(nP, **kwargs):

    dims = (1,)
    for arg in kwargs:
        # print(arg, kwargs[arg].shape)

        if isinstance(kwargs[arg], list):
            kwargs[arg] = np.array(kwargs[arg])

        if hasattr(kwargs[arg], "__len__"):
            # print(arg,kwargs[arg].shape)
            dim = kwargs[arg].shape

            """
                check if 
                    1. last dimension corresponds to number of populations
                        if not, fix!
                    2. if array has more dimensions then current dimension
                        if yes, fix!
                    3. if array is compatible with others
            """

            if not (dim[-1] == nP):
                kwargs[arg] = np.repeat(kwargs[arg][..., np.newaxis], nP, axis=-1)
                dim = kwargs[arg].shape

            if len(dim) > len(dims):
                dims = dim
            elif len(dim) == len(dims):
                for i in range(len(dim)):
                    if dim[i] > dims[i]:
                        dims = dim
                        break
        else:
            kwargs[arg] = np.atleast_1d(kwargs[arg])
            # print(arg,kwargs[arg].shape)

    return dims, kwargs
