import math
import numpy as np

class network:

    def __init__(self,**kwargs):

        para_defaults = {
            'tau_G': 0.01,
            'tau_A': 0.005,
            'tau_N': 0.2,
            'tau_M': 0.01,
            'alpha_0': 0.0,
            'J_0': -1.,

            'eta': 0.9,
            'eps': 1/np.sqrt(2),
            'alpha_0': 0.02,
            'r': 0.5,
            'kappa_E': 1.
        }

        for key in para_defaults:
            setattr(self,key,kwargs[key] if key in kwargs else para_defaults[key])

        #self.tau_G = tau_G          # synaptic timeconstant in s
        #self.tau_A = tau_A          # synaptic timeconstant in s
        #self.tau_N = tau_N          # synaptic timeconstant in s
        #self.tau_M = tau_M          # membrane timeconstant in s
        J = self.J_0 * self.tau_M          # synaptic coupling strength

        self.J = np.array([
            [self.eps * J, np.sqrt(1 - self.eps**2) * J],
            [self.eta * self.eps * J, np.sqrt(1 - (self.eta*self.eps)**2) * J]
        ])


    def I_squared_nu(self, nu, q, p):
        return - ( self.alpha(q,p)**2 + self.sigma_V(nu,p)**2 ) * np.log( (nu[p,...]/self.rate_max(nu,p))**2 * (1 + (self.alpha(q,p) / self.sigma_V(nu,p))**2) )

    def I_squared_q(self, nu, q, p):
        return -( self.alpha(q,p)**2 + 1./2 * self.sigma_V(nu, p)**2 ) * np.log( ( q[p,...]/self.rate_max(nu, p)**2 )**2 * (1 + 2*(self.alpha(q,p) / self.sigma_V(nu,p))**2) )

#    def get_q(self,nu,q,I):
#        return self.rate_max()**2 * self.sigma_V(nu) / np.sqrt(2*self.alpha(q)**2 + self.sigma_V(nu)**2) * np.exp( - I**2 / (2 * self.alpha(q)**2 + self.sigma_V(nu)**2) )

    def alpha(self, q, p=0):
        return np.sqrt(self.J[p,0]**2 * q[0,...] + self.J[p,1]**2 * q[1,...] + self.alpha_0**2)

    def sigma_V(self, nu, p=0):

        r_A = 1.-self.r
        r_N = self.r
        sigma_sq = self.J[p,0]**2 * nu[0,...] / ( 2 * (self.tau_G + self.tau_M)) + \
            (self.J[p,1]**2 * nu[1,...] / (self.tau_A + self.tau_M) * ( r_A**2/2 + (r_A * r_N * self.tau_A) / (self.tau_A + self.tau_N) ) + \
            self.J[p,1]**2 * nu[1,...] / (self.tau_N + self.tau_M) * ( r_N**2/2 + (r_A * r_N * self.tau_N) / (self.tau_A + self.tau_N) )) * self.kappa_E

        return np.sqrt(sigma_sq)

    def sigma_V_dot(self, nu, p=0):

        r_A = 1.-self.r
        r_N = self.r
        sigma_dot_sq = 1./(self.tau_G * self.tau_M) * self.J[p,0]**2 * nu[0,...] / ( 2 * (self.tau_G + self.tau_M)) + \
            (1./(self.tau_A * self.tau_M) * self.J[p,1]**2 * nu[1,...] / (self.tau_A + self.tau_M) * ( r_A**2/2 + (r_A*r_N * self.tau_A)/(self.tau_A + self.tau_N) ) + \
            1./(self.tau_N * self.tau_M) * self.J[p,1]**2 * nu[1,...] / (self.tau_N + self.tau_M) * ( r_N**2/2 + r_A * r_N * self.tau_N / (self.tau_A + self.tau_N) )) * self.kappa_E

        return np.sqrt(sigma_dot_sq)

    def rate_max(self,nu,p=0):
        return 1./(2 * math.pi) * self.sigma_V_dot(nu,p) / self.sigma_V(nu,p)

    def I(self,nu,q,p):
        return np.sqrt( np.abs(self.I_squared_nu(nu,q,p)) )

    def delta(self,nu,q,p=0):
        return self.I(nu,q,p)/self.alpha(q,p)

    def gamma(self,nu,q,p=0):
        return self.sigma_V(nu,p)/self.alpha(q,p)

    def distribution(self,nu,q,p,steps=100):

        # rate_arr = np.linspace(0,self.rate_max(),steps)
        rate_ratio = np.linspace(1/steps,1,steps-1) #rate_arr/self.rate_max()

        distr = self.gamma(nu,q,p)/(self.rate_max(nu,p)*np.sqrt(-np.pi*np.log(rate_ratio)))* \
                np.exp(-self.delta(nu,q,p)**2/2)*rate_ratio**(self.gamma(nu,q,p)**2-1)* \
                np.cosh(self.gamma(nu,q,p)*self.delta(nu,q,p)*np.sqrt(-2*np.log(rate_ratio)))

        at_zero = 0 if self.gamma(nu,q,p) > 1 else np.inf
        distr = np.insert(distr,0,at_zero)
        rate_ratio = np.linspace(0,1,steps)
        return rate_ratio, distr
