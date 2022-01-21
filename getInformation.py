import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import scipy.stats as stats
import matplotlib as mpl
from matplotlib.widgets import Slider

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *
from pythonCode.network import network

def information(steps=100,
    rateWnt=[0,20],alpha_0=[0],tau_G=[0.005],eps=[0.5],eta=[0.9],n=[0],I_alpha=[1.],I_beta=[1.],
    order=['rateWnt','alpha_0','tau_G','n','eta','eps','I_alpha','I_beta'],
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
        'I_alpha':  I_alpha,
        'I_beta':   I_beta,
        'drive':    drive,
        'mode_stats': 4,
        'J':        J,
        # 'nZeta':    51
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # return res;
    res['infoContent'][res['infoContent']==0] = np.nan
    x = np.linspace(0,1,100)
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')

    plt_para = {}
    levs = 20
    plt_para['bnw'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=levs-1)
    plt_para['heat'] = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=levs-1)
    plt_para['bnw_regions'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)


    fig = plt.figure(figsize=(8,6),dpi=300)

    ax = plt.axes([0.75,0.8,0.2,0.15])
    base  = 1.
    a = res['I_alpha'][0]
    alpha = a if a>0 else base
    beta = -a if a<0 else base
    x = np.linspace(stats.beta.ppf(0.01, alpha, beta),stats.beta.ppf(0.99, alpha, beta), 100)
    ax.plot(x,stats.beta.pdf(x,alpha,beta),label="a=%g"%a)
    plt.setp(ax,xlabel=r'$\displaystyle \nu / \nu_{max}$',
                ylabel=r'$\displaystyle \beta(x,\alpha,a)$')



    ax = plt.axes([0.025,0.1,0.6,0.75],projection="3d")
    # X = np.linspace(-3,3,options['nZeta']) ## to be replaced by data from res

    # mask_inconsistent = (regions == 3)
    # mask_no_peak = (regions == 2)
    # mask_implausible = (regions == 1)
    mask_dark_matter = (res['gamma'][0,...]**2 < 1)

    normalize = mcolors.Normalize(vmin=0, vmax=3)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['heat'])
    col_chi = s_map.to_rgba(res['chi'][0,...])

    normalize = mcolors.Normalize(vmin=0, vmax=2)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['bnw'])
    col_gamma = s_map.to_rgba(res['gamma'][0,...])

    col = col_chi
    col[mask_dark_matter,:] = col_gamma[mask_dark_matter,:]
    idx_max1 = np.nanargmax(res['infoContent'][0,...],axis=0)
    idx_max2 = np.nanargmax(res['infoContent'][0,...],axis=1)
    col[0,0,:] =[np.nan,np.nan,np.nan,np.nan]

    # col[idx_max1,range(steps),:] = [0,1,0,1]
    # col[range(steps),idx_max2,:] = [1,0,0,1]

    y_arr = res[options['order'][1]][:]
    x_arr = res[options['order'][0]][:]
    X,Y = np.meshgrid(x_arr,y_arr)
    ax.plot_surface(X,Y,res['infoContent'][0,...],facecolors=col,antialiased=False,shade=True)

    # ax.plot(res[options['order'][1]][idx_max1],res[options['order'][0]][range(steps)],res['infoContent'][0,range(steps),idx_max1]+0.1,'r',lw=3)
    # ax.plot(res[options['order'][1]][idx_max2],res[options['order'][0]][range(steps)],res['infoContent'][0,idx_max2,range(steps)]+0.001,'green',lw=3)


    plt.setp(ax,#zlim=[0,10],
        xlabel=r'$\displaystyle %s$'%options['order'][0],
        ylabel=r'$\displaystyle %s$'%options['order'][1],
        zlabel=r'$\displaystyle I_{\sum}$',
        zlim=[np.nanmin(res['infoContent']),np.nanmax(res['infoContent'])])

    # X = res[options['order'][1]][:]

    ax = plt.axes([0.75,0.45,0.2,0.2])
    # zeta_arr = np.linspace(-3,3,options['nZeta'])
    # for nu in [0.1,0.5,1,2,5]:
        # i = np.argmin(abs(nu_arr-nu))
        # print(nu,i)
    ax.plot(x_arr,res['infoContent'][0,idx_max1,range(steps)],label=r'$\displaystyle I(%s^{max})$'%(options['order'][0]))
    ax.plot(y_arr,res['infoContent'][0,range(steps),idx_max2],label=r'$\displaystyle I(%s^{max})$'%(options['order'][1]))
    ax.legend()
    plt.setp(ax,xlim=[0,0.1],
                xlabel='$\displaystyle %s \; or \; %s$'%(options['order'][0],options['order'][1]))
    # plt.setp(ax,ylim=[0,10])

    ax = plt.axes([0.75,0.1,0.2,0.2])
    # for zeta in [-3,-1.5,0,1.5,3]:
        # i = np.argmin(abs(zeta_arr-zeta))
        # print(zeta,i)
    ax.plot(x_arr,y_arr[idx_max1],label=r'$\displaystyle %s(%s^{max})$'%(options['order'][1],options['order'][0]))
    # ax.plot(x_arr[idx_max2],y_arr,label=r'$\displaystyle \%s(\%s^{max})$'%(options['order'][0],options['order'][1]))
    ax.legend()
    plt.setp(ax,xlim=[0,0.1],
                xlabel='$\displaystyle %s$'%(options['order'][0]))

    plt.suptitle('$\displaystyle \\bar{\\nu}=%.2fHz, \\alpha = %.2f, \\beta = %.2f  a=%.2f$'%(res['rateWnt'][0],res['I_alpha'][0],res['I_beta'][0],base))
    plt.show(block=False)

    # steps = 1000
    net = network()
    x = 10
    y = 10
    rate_arr , distr = net.distribution(res['rateWnt'][0],res['q'][0,x,y],steps=1000)

    return res



    fig,axes = plt.subplots(2,1,figsize=(8,12),sharex=True)
    # axes[0]
    base  = 1.
    alpha = res['I_alpha'][0]
    beta = res['I_beta'][0]
    # alpha = a if a>0 else base
    # beta = -a if a<0 else base
    # x = np.linspace(stats.beta.ppf(0.01, alpha, beta),stats.beta.ppf(0.99, alpha, beta), 100)
    axes[0].plot(rate_arr,stats.beta.pdf(rate_arr,alpha,beta),label="a=%g"%a)
    plt.setp(axes[0],xlabel=r'$\displaystyle \nu / \nu_{max}$',
                ylabel='$\displaystyle \\beta(x,\\alpha,a)$')



    p_distr, = axes[1].plot(rate_arr,distr,'k-')
    twinx = axes[1].twinx()

    p_I, = twinx.plot(rate_arr,rate_arr*distr*stats.beta.pdf(rate_arr/net.rate_max(),alpha,beta),'r-')


    def update_slider(val):
        # print('update ',val)

        val_alpha = int(samp_alpha.val/samp_alpha.valmax*steps)
        val_tau = int(samp_tau.val/samp_tau.valmax*steps)
        # print(val_alpha,val_tau)

        # val = int(val)
        _ , distr = net.distribution(res['rateWnt'][0],res['q'][0, val_alpha,val_tau],steps=1000)
        p_distr.set_ydata(distr)
        I = rate_arr*distr*stats.beta.pdf(rate_arr/net.rate_max(),alpha,beta)
        p_I.set_ydata(I)
        print('####     distribution     ####')
        print(net.rate_max())
        print(np.nanmax(distr))
        print(np.nanmax(I))
        plt.setp(axes[1],ylim=[0,np.nanmax(distr)*1.1])
        plt.setp(twinx,ylim=[0,np.nanmax(I)*1.1])
        fig.canvas.draw_idle()

    # print(steps)
    axamp_alpha = plt.axes([0.2, .03, 0.50, 0.02])
    axamp_tau = plt.axes([0.2, .06, 0.50, 0.02])

    samp_alpha = Slider(axamp_alpha, r'$\displaystyle %s$'%options['order'][0], 0, res[options['order'][0]][-1], valinit=0.0)
    samp_tau = Slider(axamp_tau, r'$\displaystyle %s$'%options['order'][1], 0, res[options['order'][1]][-1], valinit=0.005)

    update_slider(0)

    # call update function on slider value change
    samp_alpha.on_changed(update_slider)
    samp_tau.on_changed(update_slider)
    plt.setp(axes[0],xlim=[0,0.2])
    plt.setp(axes[1],xlim=[0,0.2])

    plt.show()

    return res
