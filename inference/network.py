import math
import warnings
import numpy as np
from scipy.optimize import fsolve
from .utils.structures import DistributionModelParams as distr_params


class Network:

    """
        ensure, q, sigma_V, alpha, nu_max are only calculated, not set (protected vars, attributes, ...)
    """

    def __init__(self, **kwargs):
        self.eps = kwargs.get('eps', 1/np.sqrt(2))
        self.eta = kwargs.get('eta', 0.9)

        self.populations = []
        self.n_pop = 0

    def register_population(self, idx=None,**kwargs):

        if idx is not None and (idx<self.n_pop or idx>self.n_pop):
            warnings.warn(f"Requested to add population {idx}, but currently there are {len(self.populations)} populations registered. You can only add the next population in line.", UserWarning)
            return

        self.populations.append(self.Population(len(self.populations),**kwargs))
        self.n_pop += 1

    def print_current_state(self):

        print(f"Network with {self.n_pop} populations:")
        print(f"  eps: {self.eps}")
        print(f"  eta: {self.eta}")

        print(f"  weights: {self.J}\n")
        for pop in self.populations:
            pop.print_current_state()

    def set_weights(self):
        
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
                    self.J[post.idx,pre.idx] = pre.J0 * pre.tau_M * self.eta*self.eps


    def calculate_sigma_V(self):

        for post in self.populations:

            var_V = 0.
            var_V_dot = 0.
            for pre in self.populations:
                
                for syn in pre.synapses:
                    prefactor = self.J[post.idx,pre.idx]**2 * pre.kappa * pre.nu_bar / (syn.tau_I + post.tau_M)

                    if len(pre.synapses) == 1:
                        tmp_var_V = prefactor/2.
                    elif len(pre.synapses) == 2:
                        tmp_var_V = prefactor * syn.r_I**2 / 2 + (syn.r_I * (1 - syn.r_I) * syn.tau_I) / (pre.synapses[0].tau_I + pre.synapses[1].tau_I)
                    else:
                        assert False, "not implemented for more than 2 synapses"
                    
                    var_V += tmp_var_V
                    var_V_dot += tmp_var_V / (syn.tau_I * post.tau_M)
            
            post._sigma_V = np.sqrt(var_V)
            post._nu_max = np.sqrt(var_V_dot/var_V)/(2*math.pi)


    def calculate_alpha(self):
        ## requires having set/calculated q before
        var_alpha_external = 0. ## can be implemented later
        for post in self.populations:

            var_alpha_I = 0.
            for pre in self.populations:

                var_alpha_I += self.J[post.idx,pre.idx]**2 * pre.kappa * pre.q


            post._alpha = np.sqrt(var_alpha_I + var_alpha_external + post.alpha_0**2)

    def solve_selfcon(self):
        
        q0 = np.array([pop.nu_bar**2 for pop in self.populations])
        q = fsolve(self.selfcon, q0)
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
        is_ok = all(getattr(self, attr) is not None and getattr(self, attr) > 0 for attr in attrs)
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

        def __init__(self, idx, **kwargs):

            self.idx = idx
            self.nu_bar = kwargs.get('nu_bar', 0.)
            self.alpha_0 = kwargs.get('alpha_0', 0.)
            
            self.J0 = kwargs.get('J0', -1.)

            self.kappa = kwargs.get('kappa', 1.)
            self.Psi_0 = kwargs.get('Psi_0', 0.)
            self.tau_M = kwargs.get('tau_M', 0.01)

            self.synapses = []

        @property
        def q(self):
            return self._q
        
        @q.setter
        def q(self, value):
            self._q = value

        @property
        def sigma_V(self):
            return self._sigma_V

        @property
        def alpha(self):
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


        def print_current_state(self):

            print(f"Population {self.idx}:")
            print(f"  nu_bar: {self.nu_bar}")
            print(f"  alpha_0: {self.alpha_0}")
            print(f"  q: {self.q}")
            print("")
            print(f"  alpha: {self.alpha}")
            print(f"  sigma_V: {self.sigma_V}")
            print(f"  nu_max: {self.nu_max}")
            print("")
            print(f"  tau_M: {self.tau_M}")
            print(f"  Psi_0: {self.Psi_0}")
            print(f"  kappa: {self.kappa}")
            print(f"  J0: {self.J0}")
            for syn in self.synapses:
                syn.print_current_state()

        def register_synapse(self, idx=None, **kwargs):

            if idx is not None and (idx < len(self.synapses) or idx > len(self.synapses)):
                warnings.warn(f"Requested to add synapse {idx}, but currently there are {len(self.synapses)} synapses registered. You can only add the next synapse in line.", UserWarning)
                return

            self.synapses.append(self.Synapse(len(self.synapses), **kwargs))
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

            def __init__(self, idx, **kwargs):
                self.idx = idx
                self.tau_I = kwargs.get('tau_I', 0.01)
                self.r_I = kwargs.get('r_I', 1.)
                self.tau_norm = kwargs.get('tau_norm',1.)
            
            def are_values_ok(self):
                attrs = ['tau_I','r_I','tau_norm']
                return all(getattr(self, attr) is not None and getattr(self, attr) > 0 for attr in attrs)

            def print_current_state(self):
                print(f"  Synapse {self.idx}:")
                print(f"    tau_I: {self.tau_I}")
                print(f"    r_I: {self.r_I}")
                print(f"    tau_norm: {self.tau_norm}")