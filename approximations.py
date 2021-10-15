import numpy as np
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *

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

def approximations(steps=1000,rateWnt=[0,20],alpha_0=[0,0.02,0.04],tau_G=[0.005],eps=[0.5],eta=[0.9],n=[0],J=-1.,Npop=1,drive=0,save=0,file_format='png',rerun=False,compile=False):
    
    options = {
        'Npop': Npop,
        'order': ['rateWnt','alpha_0','tau_G','n','eta','eps'],
        'rateWnt': rateWnt,
        'alpha_0': alpha_0,
        'tau_G': tau_G,
        'eps': eps,
        'eta': eta,
        'n': n,
        'drive': drive,
        'mode_stats': 2,
        'J': J
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)

    # absorb into plot_prep function
    # -----------------------------------------
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['font.size'] = 10

    steps1 = res['gamma'].shape[1]

    if len(rateWnt) > 1:
        x_key = 'rateWnt'
    if len(eps) > 1:
        x_key = 'eps'
    if len(eta) > 1:
        x_key = 'eta'
    if len(n) > 1:
        x_key = 'n'
    if len(tau_G) > 1:
        x_key = 'tau_G'

    ## define dictionary with transition point indices
    trans_idx = {}
    for key in ['inc','imp','DM','np']:
        trans_idx[key] = np.zeros(steps1,'int')
        for a in range(steps1):
            idx = np.where(res[x_key]==res[key+'_trans'][0,a])[0]
            trans_idx[key][a] = idx[0] if len(idx) else -1

    plt_para = {
        'title': {
            'descr': False,
            'x_offset': -0.2
            },
        'two_row': True
    }
    # -----------------------------------------

    lower_bound = 0.18
    box_height = 0.7
    box_width = 0.22
    v_space = 0.05
    h_space = 0.1

    fig = plt.figure(figsize=(7.5,3),dpi=300)
    ax_approx_ex = plt.axes([0.1,lower_bound,0.2,0.7])
    ax_I = plt.axes([0.2,lower_bound+0.5,0.1,0.2])
    ax_q = plt.axes([0.4,lower_bound,0.25,0.4])
    ax_q_zoom = plt.axes([0.5,lower_bound+0.5,0.15,0.2])
    ax_gamma = plt.axes([0.7,lower_bound,0.25,0.3])
    ax_chi = plt.axes([0.7,lower_bound+0.4,0.25,0.3])

    res_approx = get_approximation(0.00)
    res_approx_hetero = get_approximation(0.04)
    plot_approx_ex(ax_approx_ex)
    plot_I(ax_I,res_approx)
    plot_I(ax_I,res_approx_hetero)
    plot_q(ax_q,res,'rateWnt',trans_idx,plt_para,12,approx=True)
    plot_q_zoom(ax_q_zoom,res,'rateWnt',trans_idx,plt_para,1,approx=True)

    ax_gamma.plot(res_approx['nu'],res_approx['alpha'][:,0],'k')
    ax_gamma.plot(res_approx['nu'],res_approx['alpha'][:,1],'r')

    # plot_gamma(ax_gamma,res,'rateWnt',trans_idx,plt_para,12,approx=True)
    plot_chi(ax_chi,res,'rateWnt',trans_idx,plt_para,12,approx=True)

    plt.show(block=False)

    return res


def get_approximation(alpha=0.0):
    net = network(0.005,0.01,alpha)

    nu_lim = 12
    steps = 1001
    # I_exact = np.zeros(steps)
    res = {
        'nu': np.linspace(1/steps,nu_lim,steps),
        'I': np.zeros((steps,2)),
        'q': np.zeros((steps,2)),
        'alpha': np.zeros((steps,2)),
    }

    for i,nu in enumerate(res['nu']):
        q = np.linspace(0,max(10,4*nu**2),1001)
        min_idx = np.nanargmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))
        min_val = np.nanmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))
        # print(i, nu, min_idx, min_val)
        res['I'][i,0] = np.sqrt(net.I_squared_nu(nu,q[min_idx]))
        res['I'][i,1] = np.sqrt(net.I_squared_nu(nu,nu**2))
        res['q'][i,0] = q[min_idx]
        res['q'][i,1] = net.get_q(nu,nu**2,res['I'][i,0])
        res['alpha'][i,0] = net.alpha(res['q'][i,0])
        res['alpha'][i,1] = net.alpha(res['q'][i,1])

    return res

def plot_approx_ex(ax):

    net = network(0.005,0.01,0.00)
    ymin = -0.4
    ymax = 0
    nu = 5
    q = np.linspace(0,4*nu**2,1001)
    y_is = -np.sqrt(net.I_squared_nu(nu,nu**2)) ## intersection with function
    y_is_transformed = (y_is-ymin) / (ymax - ymin)
    # print(y_is_transformed)

    min_idx = np.nanargmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))

    ax.axvline(nu**2,ymax=y_is_transformed,color=[0.6,0.6,0.6],lw=2,ls='--')
    ax.axhline(-np.sqrt(net.I_squared_nu(nu,nu**2)),color=[0.6,0.6,0.6],lw=2,ls='--',label=r"approx.: $\displaystyle q=\bar{\nu}^2$")

    ax.plot(q,-np.sqrt(net.I_squared_nu(nu,q)),'k-',label=r"solution for $\displaystyle \bar{\nu}$")
    ax.plot(q,-np.sqrt(net.I_squared_q(nu,q)),'k--',label=r"solution for $\displaystyle q$")

    plt.setp(ax, ylim=[ymin,ymax])

def plot_I(ax,res):
    # nu_lim = 12
    # steps = res['I'].shape[0]
    ax.plot(res['nu'],(res['I'][:,0]-res['I'][:,1])/res['I'][:,0],'k',label=r"$\frac{\Delta I}{I}$")
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_exact,'k')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,0],'r')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,1],'r--')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),I_approx,'r')
