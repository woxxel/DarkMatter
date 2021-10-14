import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class network:

    def __init__(self,tau_I=0.005,tau_m=0.01,alpha_0=0.0,J=-1.):
        self.tau_I = tau_I   # synaptic timeconstant in s
        self.tau_m = tau_m   # membrane timeconstant in s
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

def plot_selfcon(nu):

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['font.size'] = 10

    q = np.linspace(0,4*nu**2,1001)

    y_range = [-0.4,0]

    net = network()

    min_idx = np.argmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))

    plt.figure(figsize=(3,2.3))
    ax = plt.axes([0.2,0.2,0.65,0.75])

    ax.plot([nu**2,nu**2],[-1,-np.sqrt(net.I_squared_nu(nu,nu**2))],color=[0.6,0.6,0.6],linewidth=3,linestyle='--')
    ax.plot(q,-np.sqrt(net.I_squared_nu(nu,np.ones(len(q))*nu**2)),'r-',linewidth=0.5,label=r"approx.: $\displaystyle q=\bar{\nu}^2$")       # plot approximation

    ax.plot(q,-np.sqrt(net.I_squared_nu(nu,q)),'k-',label=r"solution for $\displaystyle \bar{\nu}$")
    ax.plot(q,-np.sqrt(net.I_squared_q(nu,q)),'k--',label=r"solution for $\displaystyle q$")


    #plt.plot(q,net.I_squared_nu(nu,q),'k-')
    #plt.plot(q,net.I_squared_q(nu,q),'k--')
    plt.setp(ax,xticks=np.linspace(0,100,6), xticklabels=np.linspace(0,100,6).astype('int'), xlim=[0,4*nu**2])
    plt.setp(ax,yticks=np.linspace(y_range[0],y_range[1],5), yticklabels=np.linspace(y_range[0],y_range[1],5))
    ax.set_ylim(y_range)
    ax.set_xlabel(r'$\displaystyle q\,$[Hz$\displaystyle^2$]',fontsize=10)
    ax.set_ylabel(r'$\displaystyle I_0-\Psi_0$')

    ax.spines['right'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')

    ##print q[min_idx]
    ##print np.sqrt(net.I_squared_q(nu,q[min_idx]))
    ax.annotate(r'$\displaystyle q = \bar{\nu}^2$',xy=[nu**2,-np.sqrt(net.I_squared_nu(nu,nu**2))],xytext=[nu**2-22,-0.15],arrowprops=dict(arrowstyle="->"),fontsize=10)
    ax.annotate(r'$\displaystyle (\bar{\nu}^{\star},q^{\star})$',xy=[q[min_idx],-np.sqrt(net.I_squared_nu(nu,q[min_idx]))],xytext=[nu**2+15,-0.3],arrowprops=dict(arrowstyle="->"),fontsize=10)

    ax.legend(prop={'size':12},frameon=False,loc=[0.37,0.65])

    # sv_name = '../paper draft/inhib only/pics/selfcon_approx.png'
    # print('Figure saved as "%s"' % sv_name)
    # plt.savefig(sv_name,dpi=300)

    plt.show(block=False)
