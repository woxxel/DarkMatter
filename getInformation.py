import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import math
import scipy.stats as stats
import matplotlib as mpl

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *

def information(steps=100,
    rateWnt=[0,20],alpha_0=[0],tau_G=[0.005],eps=[0.5],eta=[0.9],n=[0],zeta=[-3],
    order=['rateWnt','alpha_0','tau_G','n','eta','eps','zeta'],
    J=-1.,Npop=1,drive=0,
    save=0,file_format='png',
    rerun=False,compile=False):

    options = {
        'Npop':     Npop,
        'order':    order,
        'rateWnt':  rateWnt,
        'alpha_0':  alpha_0,
        'tau_G':    tau_G,
        'eps':      eps,
        'eta':      eta,
        'n':        n,
        'zeta':     zeta,
        'drive':    drive,
        'mode_stats': 4,
        'J':        J,
        # 'nZeta':    51
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # return res;
    x = np.linspace(0,1,100)
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')

    fig = plt.figure(figsize=(8,6),dpi=300)

    ax = plt.axes([0.1,0.65,0.3,0.25])
    for a in np.linspace(-3,3,5):
        alpha = a if a>0 else 1
        beta = -a if a<0 else 1
        ax.plot(x,stats.beta.pdf(x,alpha,beta),label="a=%g"%a)
    ax.legend();plt.show(block=False)
    ax.set_xlabel(r'$\displaystyle \nu / \nu_{max}$')


    ax = plt.axes([0.5,0.6,0.45,0.35],projection="3d")
    # X = np.linspace(-3,3,options['nZeta']) ## to be replaced by data from res
    Y = res[options['order'][0]][:]
    X = res[options['order'][1]][:]
    X,Y = np.meshgrid(X,Y)
    ax.plot_surface(X,Y,res['infoContent'][0,...])
    plt.setp(ax,#zlim=[0,10],
        xlabel=r'$\displaystyle %s$'%options['order'][1],
        ylabel=r'$\displaystyle %s$'%options['order'][0],
        zlabel=r'$\displaystyle I_{\sum}$')

    ax = plt.axes([0.6,0.1,0.35,0.4])
    nu_arr = res[options['order'][0]][:]
    zeta_arr = res[options['order'][1]][:]
    # zeta_arr = np.linspace(-3,3,options['nZeta'])
    for nu in [0.1,0.5,1,2,5]:
        i = np.argmin(abs(nu_arr-nu))
        print(nu,i)
        ax.plot(zeta_arr,res['infoContent'][0,i,:],label=r'$\displaystyle \nu=%.2g$'%nu)
    ax.legend()
    # plt.setp(ax,ylim=[0,10])

    ax = plt.axes([0.1,0.1,0.35,0.4])
    for zeta in [-3,-1.5,0,1.5,3]:
        i = np.argmin(abs(zeta_arr-zeta))
        print(zeta,i)
        ax.plot(nu_arr,res['infoContent'][0,:,i],label=r'$\displaystyle a=%.2g$'%zeta)
    ax.legend()
    # plt.setp(ax,ylim=[0,10])

    plt.show(block=False)



    return res
