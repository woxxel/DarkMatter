import numpy as np
import matplotlib.pyplot as plt


def spike_amplitude(steps=100):
    
    n_arr = np.linspace(0,1,steps)
        
    x = np.zeros(steps)
    y = np.zeros(steps)
    
    tau_A = 0.005
    tau_N = 0.2
    
    for i in range(steps):
        n = n_arr[i]
        a = 1-n_arr[i]
        
        x[i] = n
        y[i] = a*1./tau_A + n*1./tau_N
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.figure()
    plt.plot(x,y)
    plt.xlabel('n')
    plt.ylabel('amplitude')
    plt.grid()
    plt.show(block=False)

class tempVar:
    
    tau_A = 0.005
    tau_N = 0.2
    tau_M = 0.01
    
    eps = np.sqrt(0.5)
    eta = 0.9
    
    
    nu_I = 1
    nu_E = 1
    
    tau_G = 0.03
    n = 0.2
    a = 1-n
    
    def __init__(self,nu,eps,eta):
        self.nu_I = nu
        self.nu_E = nu
        
        self.J_IE = eps**2*self.tau_M
        self.J_II = np.sqrt(1-eps**2)**2*self.tau_M
    
    
    def sigma_E(self):
        
        sigma_tmp = self.J_IE**2*self.nu_E * (1./(self.tau_A + self.tau_M) * (self.a**2/2. + self.a*self.n*self.tau_A/(self.tau_A + self.tau_N)) + 1./(self.tau_N + self.tau_M) * (self.n**2/2. + self.a*self.n*self.tau_N/(self.tau_A + self.tau_N)))
        return np.sqrt(sigma_tmp)
        #print stringPrint
        
    def sigma_I(self):
    
        sigma_tmp = self.J_II**2*self.nu_I * (1./(2*(self.tau_G + self.tau_M)))
        return np.sqrt(sigma_tmp)
        
    def print_var(self):
        
        print "Variance sigma_E: ", self.sigma_E(), ",\t Variance sigma_I: ", self.sigma_I()
        
    def n_loop(self,steps=100):
        n = np.linspace(0,1,steps)
        
        for i in range(steps):
            self.n = n[i]
            self.a = 1-n[i]
            self.print_var()
        
    def find_n(self,steps):
        n = np.linspace(0,1,steps)
        
        x = np.zeros(steps)
        y = np.zeros(steps)
        
        for i in range(steps):
            self.n = n[i]
            self.a = 1-n[i]
            tau_G_tmp = (2*self.sigma_E()**2/(self.nu_I*self.J_II**2))**(-1) - self.tau_M
            
            x[i] = n[i]
            y[i] = tau_G_tmp
            print "n = %4.2f, \t tau_G = %g" %(self.n,tau_G_tmp)
        
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
  
        plt.figure()
        plt.plot(x,y)
        plt.xlabel('n')
        plt.ylabel(r'$\displaystyle \tau_G [ms]$')
        plt.grid()
        plt.show(block=False)
        
#def get_tau_G_n():
    #a = tempVar