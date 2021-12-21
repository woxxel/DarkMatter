import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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

    plt_para = {}
    levs = 20
    plt_para['bnw'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=levs-1)
    plt_para['heat'] = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=levs-1)
    plt_para['bnw_regions'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)


    fig = plt.figure(figsize=(8,6),dpi=300)

    ax = plt.axes([0.1,0.65,0.3,0.25])
    base  = 1.
    # for a in np.linspace(-3,3,7):
    #     alpha = a if a>0 else base
    #     beta = -a if a<0 else base
    #     x = np.linspace(stats.beta.ppf(0.01, alpha, beta),stats.beta.ppf(0.99, alpha, beta), 100)
    #     ax.plot(x,stats.beta.pdf(x,alpha,beta),label="a=%g"%a)
    # ax.legend();plt.show(block=False)
    a = res['zeta'][0]
    alpha = a if a>0 else base
    beta = -a if a<0 else base
    x = np.linspace(stats.beta.ppf(0.01, alpha, beta),stats.beta.ppf(0.99, alpha, beta), 100)
    ax.plot(x,stats.beta.pdf(x,alpha,beta),label="a=%g"%a)
    ax.set_xlabel(r'$\displaystyle \nu / \nu_{max}$')



    ax = plt.axes([0.5,0.6,0.45,0.35],projection="3d")
    # X = np.linspace(-3,3,options['nZeta']) ## to be replaced by data from res
    normalize = mcolors.Normalize(vmin=0, vmax=3)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['heat'])
    col_chi = s_map.to_rgba(res['chi'][0,...])

    normalize = mcolors.Normalize(vmin=0, vmax=2)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['bnw'])
    col_gamma = s_map.to_rgba(res['gamma'][0,...])
    print(col_chi.shape)
    Y = res[options['order'][0]][:]
    X = res[options['order'][1]][:]
    X,Y = np.meshgrid(X,Y)
    # ax.plot_surface(X,Y,res['infoContent'][0,...],facecolors=col_gamma,antialiased=False,shade=False)
    ax.plot_surface(X,Y,res['infoContent'][0,...],facecolors=col_chi,antialiased=False,shade=False)
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
