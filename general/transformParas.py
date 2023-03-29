import math
import numpy as np
from .helper import *

class transformParas:

    def __init__(self,gamma=1.5,delta=4.2,nu_max=45.,
                 J=-1,tau_M=0.01,tau_A=0.005,tau_N=0.2,tau_G=0.005,
                 nu_E=None,nu_I=2.,kappa_E=4.):

        self.gamma = gamma
        self.delta = delta
        self.nu_max = nu_max

        t_scale = 1.
        self.tau_M = tau_M/t_scale
        self.tau_A = tau_A/t_scale
        self.tau_N = tau_N/t_scale
        self.tau_G = tau_G/t_scale

        self.J_0 = J * tau_M

        self.nu_E = nu_E
        self.nu_I = nu_I

        self.kappa_E = kappa_E

        self.calculated = np.zeros(3,dtype='bool')


    def updateParas(self,newParas):

        for key,val in newParas.items():
            if hasattr(self,key):
                setattr(self,key,val)

        self.calculated[:] = False


    def get_nu_mean(self):
        if self.nu_E:
            self.nu_mean = self.nu_E
            return self.nu_E
        elif self.calculated[0]:
            return self.nu_mean
        else:
            self.nu_mean = self.nu_max * self.gamma / np.sqrt(self.gamma**2 + 1) * np.exp(-self.delta**2 / (2*(1+self.gamma**2)))
            self.calculated[0] = True
            return self.nu_mean

    def get_alpha_0(self):
        if self.calculated[1]:
            return self.alpha_0
        else:

            nu_mean = self.get_nu_mean()

            q = self.nu_max**2 * self.gamma / np.sqrt(self.gamma**2 + 2) * np.exp(-self.delta**2 / (2+self.gamma**2))
            self.alpha_0 = np.sqrt(self.J_0**2 * ( self.nu_mean/ (2 * self.gamma**2 * (self.tau_A + self.tau_M)) - q))
            self.calculated[1] = True
            return self.alpha_0

    def get_r_NMDA(self):
        if self.calculated[2]:
            return self.r_NMDA
        else:
            rFinder = RootFinder(-1, 3, 0.01)
            roots = rFinder.find(lambda r: self.get_nu_max(r) - self.nu_max)
            roots = [ro for ro in roots if 0<=ro<=1]

            self.r_NMDA = roots[0] if len(roots) else np.NaN
            self.calculated[2] = True
            return self.r_NMDA

    def returnParas(self,paras={}):
        self.updateParas({'gamma':paras[0],'delta':paras[1],'nu_max':paras[2]})

        if not self.calculated[0]:
            self.get_nu_mean()
        if not self.calculated[1]:
            self.get_alpha_0()
        if not self.calculated[2]:
            self.get_r_NMDA()

        return self.nu_mean, self.alpha_0, self.r_NMDA



    def sigma_V(self, r = 0.5):

        r_A = 1.-r
        r_N = r

        nu_E = self.get_nu_mean()
        #nu_I = nu_E*self.nu_I

        sigma_V_G = self.J_0**2 * self.nu_I / ( 2 * (self.tau_G + self.tau_M))

        sigma_V_A = self.J_0**2 * nu_E / (self.tau_A + self.tau_M) * \
         ( r_A**2/2. + (r_A * r_N * self.tau_A) / (self.tau_A + self.tau_N) )

        sigma_V_N = self.J_0**2 * nu_E / (self.tau_N + self.tau_M) * \
         ( r_N**2/2. + (r_A * r_N * self.tau_N) / (self.tau_A + self.tau_N) )

        return np.sqrt(sigma_V_G + (sigma_V_A + sigma_V_N) * self.kappa_E)

    def sigma_V_dot(self, r = 0.5):

        r_A = 1.-r
        r_N = r

        nu_E = self.get_nu_mean()
        #nu_I = nu_E*self.nu_I

        sigma_dV_G = 1./(self.tau_G * self.tau_M) * self.J_0**2 * self.nu_I / ( 2 * (self.tau_G + self.tau_M))

        sigma_dV_A = 1./(self.tau_A * self.tau_M) * self.J_0**2 * nu_E / (self.tau_A + self.tau_M) * \
         ( r_A**2/2. + (r_A * r_N * self.tau_A) / (self.tau_A + self.tau_N) )

        sigma_dV_N = 1./(self.tau_N * self.tau_M) * self.J_0**2 * nu_E / (self.tau_N + self.tau_M) * \
         ( r_N**2/2. + (r_A * r_N * self.tau_N) / (self.tau_A + self.tau_N) )

        return np.sqrt(sigma_dV_G + (sigma_dV_A + sigma_dV_N) * self.kappa_E)


    def get_nu_max(self,r):
        res = 1./(2*math.pi) * self.sigma_V_dot(r)/self.sigma_V(r)
        #print('paras:',r,self.sigma_V(r),self.sigma_V_dot(r))
        #print(r,res)
        return res
