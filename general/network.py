import math
import warnings
import numpy as np
from scipy.optimize import fsolve
from .utils.structures import DistributionModelParams as distr_params

# import inspect


class Network:

    """
        ensure, q, sigma_V, alpha, nu_max are only calculated, not set (protected vars, attributes, ...)
    """
    n_pop = 0

    _recalculate_sigma_V = True
    _recalculate_alpha = True

    _counter_sigma_V = 0
    _counter_alpha = 0

    def __init__(self, **kwargs):
        self.populations = []

        self._eps = kwargs.get('eps', 1/np.sqrt(2))
        self._eta = kwargs.get('eta', 0.9)

    @property
    def eps(self):
        return self._eps

    @eps.setter
    def eps(self, value):
        self._eps = value
        self._set_weights()

    @property
    def eta(self):
        return self._eta

    @eta.setter
    def eta(self, value):
        self._eta = value
        self._set_weights()

    def set_param(self,key,var):
        name, indices = parse_name_and_indices(key, ["p","s"])
        p = indices[0]
        s = indices[1]
        if p is None:
            assert hasattr(self, name), f"Network has no attribute {name}!"
            setattr(self, name, var)
        else:
            pop = self.populations[p]
            if s is None:
                assert hasattr(pop, name), f"Population has no attribute {name}!"
                setattr(pop, name, var)
            else:
                syn = pop.synapses[s]
                assert hasattr(syn, name), f"Synapse has no attribute {name}!"
                setattr(syn, name, var)
            

    def register_population(self, idx=None,**kwargs):

        if idx is not None and (idx<self.n_pop or idx>self.n_pop):
            warnings.warn(f"Requested to add population {idx}, but currently there are {len(self.populations)} populations registered. You can only add the next population in line.", UserWarning)
            return
        ## check if idx==len(self.populations)

        self.populations.append(Population(self,**kwargs))
        self.n_pop += 1

        self._set_weights()

    def print_current_state(self):

        print(f"Network with {self.n_pop} populations:")
        print(f"  eps: {self.eps}")
        print(f"  eta: {self.eta}")

        print(f"  weights: {self.J}\n")
        for pop in self.populations:
            pop.print_current_state()

    def _set_weights(self):

        self.J = np.zeros((self.n_pop,self.n_pop))
        for post in self.populations:
            for pre in self.populations:
                if post.J0 < 0 and pre.J0 < 0:      # J_II
                    self.J[post.idx,pre.idx] = pre.J0 * pre.tau_M * np.sqrt(1 - self.eps**2)
                elif post.J0 > 0 and pre.J0 < 0:    # J_EI
                    self.J[post.idx,pre.idx] = pre.J0 * pre.tau_M * np.sqrt(1 - (self.eta*self.eps)**2)
                elif post.J0 < 0 and pre.J0 > 0:    # J_IE
                    self.J[post.idx,pre.idx] = pre.J0 * pre.tau_M * self.eps
                elif post.J0 > 0 and pre.J0 > 0:    # J_EE
                    self.J[post.idx,pre.idx] = pre.J0 * pre.tau_M * self.eta * self.eps
                # print(f"J(p={post.idx},{pre.idx})={self.J[post.idx,pre.idx]}")

        self._recalculate_sigma_V = True
        self._recalculate_alpha = True

    def calculate_sigma_V(self):

        if not self._recalculate_sigma_V:
            return

        for post in self.populations:
            var_V = 0.
            var_V_dot = 0.
            for pre in self.populations:
                # print(f"in post/pre {post.idx},{pre.idx}")

                for syn in pre.synapses:
                    prefactor = self.J[post.idx,pre.idx]**2 * pre.kappa * pre.nu_bar / (syn.tau_I + post.tau_M)
                    # print(f"J: {self.J[post.idx,pre.idx]}, tau_I: {syn.tau_I}, tau_M: {post.tau_M}")
                    # print(f"prefactor: {prefactor}")
                    if len(pre.synapses) == 1:
                        tmp_var_V = prefactor/2.
                    elif len(pre.synapses) == 2:
                        tmp_var_V = prefactor * (syn.r_I**2 / 2 + (syn.r_I * (1 - syn.r_I) * syn.tau_I) / (pre.synapses[0].tau_I + pre.synapses[1].tau_I))
                    else:
                        assert False, "not implemented for more than 2 synapses"
                    # print(f"tmp_var_V: {tmp_var_V}")
                    var_V += tmp_var_V
                    var_V_dot += tmp_var_V / (syn.tau_I * post.tau_M)

            post._sigma_V = np.sqrt(var_V)
            post._nu_max = (
                np.sqrt(var_V_dot / var_V) / (2 * math.pi) if var_V > 0 else np.nan
            )
        self._counter_sigma_V += 1
        self._recalculate_sigma_V = False

    def calculate_alpha(self):

        if not self._recalculate_alpha:
            return

        ## requires having set/calculated q before
        var_alpha_external = 0. ## can be implemented later
        for post in self.populations:

            var_alpha_I = 0.
            for pre in self.populations:
                # print(f"in post/pre {post.idx},{pre.idx}")
                # print(f"J: {self.J[post.idx,pre.idx]}, kappa: {pre.kappa}, q: {pre.q}, alpha_0: {post.alpha_0}")
                var_alpha_I += self.J[post.idx,pre.idx]**2 * pre.kappa * pre.q
                # print(f"var_alpha_I: {var_alpha_I}")

            post._alpha = np.sqrt(var_alpha_I + var_alpha_external + post.alpha_0**2)

        self._counter_alpha += 1
        self._recalculate_alpha = False

    def solve_selfcon(self):

        q0 = np.array([pop.nu_bar**2 for pop in self.populations])
        q = fsolve(self.selfcon, q0, xtol=1e-12)
        for pop in self.populations:
            pop.q = q[pop.idx]

        self.calculate_alpha()

        return np.allclose(self.selfcon(q),0)#, "selfconsistency not achieved!"

    def selfcon(self,q):

        for pop in self.populations:
            pop.q = q[pop.idx]

        self.calculate_alpha()

        dI = np.zeros(self.n_pop)
        for pop in self.populations:
            dI[pop.idx] = pop.selfcon()

        return dI

    def are_values_ok(self):

        attrs = ["eps", "eta"]
        is_ok = all(getattr(self, attr) is not None and getattr(self, attr) >= 0 for attr in attrs)
        if not is_ok:
            return False

        for pop in self.populations:
            if not pop.are_values_ok():
                return False
        return True

    def export_distr(self):

        distr = []
        for pop in self.populations:
            distr.append(distr_params(pop.gamma,pop.delta,pop.nu_max,pop.kappa))
        return distr

class Population:

    _sigma_V = 0.
    _alpha = 0.
    _nu_max = 0.
    _q = 0.
    

    def __init__(self, network, **kwargs):

        # print("IMPLEMENT CALCULATING (approximating?) IMPLAUSIBLE REGION")
        self.synapses = []

        # assert isinstance(network, Network), "Population has to be part of a Network class instance!"
        self.network = network
        self.idx = len(network.populations)

        self.nu_bar = kwargs.get('nu_bar', 0.)
        self.alpha_0 = kwargs.get('alpha_0', 0.)
        
        self._J0 = kwargs.get('J0', -1.)

        self.kappa = kwargs.get('kappa', 1.)
        self.Psi_0 = kwargs.get('Psi_0', 0.)
        self._tau_M = kwargs.get('tau_M', 0.01)

    @property
    def nu_bar(self):
        return self._nu_bar

    @nu_bar.setter
    def nu_bar(self, value):
        self._nu_bar = value
        self.network._recalculate_sigma_V = True

    @property
    def alpha_0(self):
        return self._alpha_0
    
    @alpha_0.setter
    def alpha_0(self, value):
        self._alpha_0 = value
        self.network._recalculate_alpha = True

    @property
    def tau_M(self):
        return self._tau_M
    
    @tau_M.setter
    def tau_M(self, value):
        self._tau_M = value
        self.network._set_weights()

    @property
    def J0(self):
        return self._J0

    @J0.setter
    def J0(self, value):
        self._J0 = value
        self.network._set_weights()

    @property
    def q(self):
        return self._q
    
    @q.setter
    def q(self, value):
        self._q = value
        self.network._recalculate_alpha = True

    @property
    def sigma_V(self):
        self.network.calculate_sigma_V()
        return self._sigma_V

    @property
    def alpha(self):
        self.network.calculate_alpha()
        return self._alpha

    @property
    def nu_max(self):
        return self._nu_max
    
    @property
    def gamma(self):
        return self.sigma_V / self.alpha
    
    @property
    def delta(self):
        return self.I_nu_squared()**0.5 / self.alpha

    @property
    def chi(self):
        return -np.log10(self.nu_peak()/self.nu_bar)
    
    
    def are_values_ok(self):
        attrs = ['nu_bar','alpha_0','kappa','tau_M']
        is_ok = all(getattr(self, attr) is not None and getattr(self, attr) > 0 for attr in attrs)
        if not is_ok:
            return False

        for syn in self.synapses:
            if not syn.are_values_ok():
                return False
        return True

    def distribution(self,steps=1000):

        # rate_arr = np.linspace(0,self.nu_max(),steps)
        rate_ratio = np.linspace(1/steps,1,steps-1) #rate_arr/self.nu_max()

        distr = self.gamma/(self.nu_max*np.sqrt(-np.pi*np.log(rate_ratio)))* \
                np.exp(-self.delta**2/2)*rate_ratio**(self.gamma**2-1)* \
                np.cosh(self.gamma*self.delta*np.sqrt(-2*np.log(rate_ratio)))

        at_zero = 0 if self.gamma > 1 else np.inf
        distr = np.insert(distr,0,at_zero)
        rate_ratio = np.linspace(0,1,steps)
        return rate_ratio, distr


    def nu_peak(self):
        return self.nu_max * np.exp( - (self.gamma**2 * self.delta**2 - 2*(self.gamma**2 - 1) + self.gamma * self.delta *np.sqrt(self.gamma**2 * self.delta**2 - 4*(self.gamma**2 - 1))) / (4 * (self.gamma**2 - 1)**2))


    def is_selfconsistent_nu(self):
        return (self.nu_bar/self.nu_max)**2 * (self.alpha**2 + self.sigma_V**2)/self.sigma_V**2 < 1
    
    def is_selfconsistent_q(self):
        return (self.q/self.nu_max**2)**2 * (2 * self.alpha**2 + self.sigma_V**2)/self.sigma_V**2 < 1
    
    @property
    def is_selfconsistent(self):
        return self.is_selfconsistent_nu() and self.is_selfconsistent_q()

    @property
    def is_dark_matter(self):
        return self.gamma > 1

    @property
    def is_no_peak(self):
        return (self.gamma * self.delta)**2 - 4*(self.gamma**2 - 1) < 0

    @property
    def is_implausible(self):
        ## this is an approximation - implement integration?
        n_steps = 1000
        rho = self.distribution(n_steps)[1]
        rho = rho[np.isfinite(rho)]
        n_steps = len(rho)
        return rho[int(0.9 * n_steps):].sum()/rho.sum() > 0.1

    def print_current_state(self):

        print(f"Population {self.idx}:")
        print(f"  nu_bar: {self.nu_bar}")
        print(f"  alpha_0: {self.alpha_0}")
        print(f"  q: {self.q}")
        print("")
        print(f"  alpha: {self.alpha}")
        print(f"  sigma_V: {self.sigma_V}")
        print("")
        print(f"  tau_M: {self.tau_M}")
        print(f"  Psi_0: {self.Psi_0}")
        print(f"  kappa: {self.kappa}")
        print(f"  J0: {self.J0}")
        for syn in self.synapses:
            syn.print_current_state()

    def print_meta_parameter(self):
        print(f"  Meta parameters for population {self.idx}:")
        print(f"    gamma: {self.gamma}")
        print(f"    delta: {self.delta}")
        print(f"    nu_max: {self.nu_max}")

    def register_synapse(self, idx=None, **kwargs):

        if idx is not None and (idx < len(self.synapses) or idx > len(self.synapses)):
            warnings.warn(f"Requested to add synapse {idx}, but currently there are {len(self.synapses)} synapses registered. You can only add the next synapse in line.", UserWarning)
            return

        self.synapses.append(Synapse(self, **kwargs))
        self.n_synapses = len(self.synapses)

        # self.synapses.append(self.Synapse(len(self.synapses), **kwargs))
        # pass
    

    def check_mixture(self):

        r = 0
        for syn in self.synapses:
            r += syn.r_I
        
        assert r==1, "synapse ratio within population {p} has to add up to 1!"
    

    def I_nu_squared(self):

        return - (self.alpha**2 + self.sigma_V**2) * np.log( (self.nu_bar/self.nu_max)**2 * (1 + (self.alpha/self.sigma_V)**2) )

    
    def I_q_squared(self):

        return - (self.alpha**2 + 1./2 * self.sigma_V**2) * np.log( (self.q/self.nu_max**2)**2 * (1 + 2*(self.alpha/self.sigma_V)**2) )

    
    def selfcon(self):
        return self.I_nu_squared() - self.I_q_squared() 


class Synapse:

    def __init__(self, population, **kwargs):

        self.population = population
        self.idx = len(population.synapses)
        self._tau_I = kwargs.get('tau_I', 0.01)
        self._r_I = kwargs.get('r_I', 1.)
        self._tau_norm = kwargs.get('tau_norm',1.)

    @property
    def tau_I(self):
        return self._tau_I
    
    @tau_I.setter
    def tau_I(self, value):
        self._tau_I = value
        self.population.network._recalculate_sigma_V = True

    @property
    def r_I(self):
        return self._r_I

    @r_I.setter
    def r_I(self, value):
        self._r_I = value
        self.population.network._recalculate_sigma_V = True

    @property
    def tau_norm(self):
        return self._tau_norm

    @tau_norm.setter
    def tau_norm(self, value):
        self._tau_norm = value
        self.population.network._recalculate_sigma_V = True

    def are_values_ok(self):
        attrs = ['tau_I','r_I','tau_norm']
        return all(getattr(self, attr) is not None and getattr(self, attr) > 0 for attr in attrs)

    def print_current_state(self):
        print(f"  Synapse {self.idx}:")
        print(f"    tau_I: {self.tau_I}")
        print(f"    r_I: {self.r_I}")
        print(f"    tau_norm: {self.tau_norm}")



def get_nu_bar(param: distr_params):
    return param.nu_max * param.gamma / np.sqrt(param.gamma**2 + 1) * np.exp(-param.delta**2 / (2*(1+param.gamma**2)))

def get_q(param: distr_params):
    return param.nu_max**2 * param.gamma / np.sqrt(param.gamma**2 + 2) * np.exp(-param.delta**2 / (2+param.gamma**2))


import re
from typing import Iterable, List, Optional, Tuple

def parse_name_and_indices(
    s: str, literals: Iterable[str]
) -> Tuple[str, List[Optional[int]]]:
    """
    Returns (variable_name, [idx_or_None per literal in the same order]).
    Variable name = prefix before the first <literal><digits> token.
    """

    if isinstance(literals, str):
        lits = literals.split()
    else:
        lits = list(literals)

    alts = "|".join(map(re.escape, lits))
    # Don't match inside letter-words; allow underscores and punctuation as separators.
    rx = re.compile(rf"_(?<![A-Za-z])({alts})(\d+)(?![A-Za-z])")

    found = {}
    first_pos = None

    for m in rx.finditer(s):
        lit, num = m.group(1), int(m.group(2))
        if lit not in found:  # keep only first per literal
            found[lit] = (m.start(), num)
            if first_pos is None or m.start() < first_pos:
                first_pos = m.start()

    # Remove all matches from the string and output the remaining string
    name = rx.sub("", s)

    # name = s[: first_pos - 1] if first_pos is not None else s
    indices: List[Optional[int]] = [
        found[lit][1] if lit in found else None for lit in lits
    ]
    return name, indices