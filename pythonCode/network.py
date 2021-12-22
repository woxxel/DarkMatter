import math
import numpy as np

class network:

    def __init__(self,tau_I=0.005,tau_m=0.01,alpha_0=0.0,J=-1.):
        self.tau_I = tau_I          # synaptic timeconstant in s
        self.tau_m = tau_m          # membrane timeconstant in s
        self.J = J * tau_m          # synaptic coupling strength
        self.alpha_0 = alpha_0

    def I_squared_nu(self, nu, q):
        return - ( self.alpha(q)**2 + self.sigma_V(nu)**2 ) * np.log( (nu/self.rate_max())**2 * (1 + (self.alpha(q) / self.sigma_V(nu))**2) )

    def I_squared_q(self, nu, q):
        return -( self.alpha(q)**2 + 1./2 * self.sigma_V(nu)**2 ) * np.log( ( q/self.rate_max()**2 )**2 * (1 + 2*(self.alpha(q) / self.sigma_V(nu))**2) )

    def get_q(self,nu,q,I):
        return self.rate_max()**2 * self.sigma_V(nu) / np.sqrt(2*self.alpha(q)**2 + self.sigma_V(nu)**2) * np.exp( - I**2 / (2 * self.alpha(q)**2 + self.sigma_V(nu)**2) )

    def alpha(self, q):
        return np.sqrt(self.J**2 * q + self.alpha_0**2)

    def sigma_V(self, nu):
        return np.sqrt((self.J**2 * nu) / self.tau_q())

    def rate_max(self):
        return (2 * math.pi * np.sqrt(self.tau_I*self.tau_m))**(-1)

    def tau_q(self):
        return 2 * (self.tau_I + self.tau_m)

    def I(self,nu,q):
        return np.sqrt( self.I_squared_nu(nu,q) )

    def delta(self,nu,q):
        return self.I(nu,q)/self.alpha(q)

    def gamma(self,nu,q):
        return self.sigma_V(nu)/self.alpha(q)

    def distribution(self,nu,q,steps=100):

        # rate_arr = np.linspace(0,self.rate_max(),steps)
        rate_ratio = np.linspace(0,1,steps)#rate_arr/self.rate_max()

        distr = self.gamma(nu,q)/(self.rate_max()*np.sqrt(-np.pi*np.log(rate_ratio)))*np.exp(-self.delta(nu,q)**2/2)*rate_ratio**(self.gamma(nu,q)**2-1)*np.cosh(self.gamma(nu,q)*self.delta(nu,q)*np.sqrt(-2*np.log(rate_ratio)))

        return rate_ratio, distr
