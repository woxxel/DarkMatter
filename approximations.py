import numpy as np
from numpy.ma import masked_array
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *
from pythonCode.network import network

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
    # res = darkMatter(steps=steps,options=options)


    options = {
        'Npop': Npop,
        'order': ['rateWnt','alpha_0','tau_G','n','eta','eps'],
        'rateWnt': [0,20],
        'alpha_0': [0,0.2],
        'tau_G': tau_G,
        'eps': eps,
        'eta': eta,
        'n': n,
        'drive': drive,
        'mode_stats': 3,
        'J': J
    }
    res_sharks = darkMatter(steps=500,options=options,rerun=rerun,compile=compile)


    # absorb into plot_prep function
    # -----------------------------------------
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')

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

    upper_row = 0.55
    lower_row = 0.1
    box_height = 0.75
    box_width = 0.22
    v_space = 0.05
    h_space = 0.1

    fig = plt.figure(figsize=(7.5,5),dpi=300)
    ax_approx_ex = plt.axes([0.1,upper_row,0.225,0.39])
    ax_I = plt.axes([0.25,upper_row+0.325,0.125,0.1])
    ax_q = plt.axes([0.45,upper_row,0.225,0.2])
    ax_q_zoom = plt.axes([0.55,upper_row+0.325,0.125,0.1])
    ax_alpha = plt.axes([0.775,upper_row+0.28,0.2,0.115])
    ax_gamma = plt.axes([0.775,upper_row+0.14,0.2,0.115])
    ax_chi = plt.axes([0.775,upper_row,0.2,0.115])

    ax_sharkfin = plt.axes([0.1,lower_row,0.2,0.3])
    ax_sharkfin_approx = plt.axes([0.325,lower_row,0.2,0.3])
    ax_sharkfin_KL = plt.axes([0.675,lower_row,0.2,0.3])

    plot_approx(ax_approx_ex,ax_I,ax_alpha,res)
    set_title(ax_approx_ex,1,'',(-0.15,0),10)

    # set_title(ax_approx_ex,order=1,title=r'approximation of $\displaystyle I$')

    plot_q(ax_q,res,'rateWnt',trans_idx,plt_para,12,approx=True,order=0)
    set_title(ax_q,2,'',(-0.15,0),70)
    ax_q.plot(np.NaN,np.NaN,color='gray',ls='--',label='approx.')
    ax_q.legend(prop={'size':10},handlelength=1,bbox_to_anchor=(0.0,1.2),loc='upper left',frameon=False)
    plot_q_zoom(ax_q_zoom,res,'rateWnt',trans_idx,plt_para,1,idxs=[0,2],approx=True,order=0)

    plot_gamma(ax_gamma,res,'rateWnt',trans_idx,plt_para,12,idxs=[0,2],approx=True,order=0)
    plot_chi(ax_chi,res,'rateWnt',trans_idx,plt_para,12,idxs=[0,2],approx=True,order=0)
    set_title(ax_alpha,3,'',(-0.15,0),-20)

    for ax in [ax_approx_ex,ax_I,ax_q,ax_q_zoom,ax_alpha,ax_gamma,ax_chi]:
        remove_frame(ax,['top','right'])

    plt.setp(ax_approx_ex,ylabel=r'$\displaystyle \Psi_0 - I_0$',xlabel=r'$\displaystyle q$ [Hz$^2$]')
    plt.setp(ax_alpha,xticklabels=[],xlabel=None,ylabel=r'$\displaystyle \alpha$')
    plt.setp(ax_gamma,xticklabels=[],xlabel=None,ylabel=r'$\displaystyle \gamma$')
    plt.setp(ax_chi,ylabel=r'$\displaystyle \chi$')

    p = 0
    _ = plot_fins(ax_sharkfin,res_sharks[options['order'][0]],res_sharks[options['order'][1]],res_sharks['gamma'][p,...],res_sharks['chi'][p,...],res_sharks['regions'][p,...],plt_para)
    plt.setp(ax_sharkfin,ylabel=r'$\displaystyle \alpha_0$',xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]')
    pchi,pgamma = plot_fins(ax_sharkfin_approx,res_sharks[options['order'][0]],res_sharks[options['order'][1]],res_sharks['gamma_approx'][p,...],res_sharks['chi_approx'][p,...],res_sharks['regions_approx'][p,...],plt_para)
    plt.setp(ax_sharkfin_approx,yticklabels=[],xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]')
    plot_colorbar(pchi,pgamma,plt_para,x=[0.55,0.57],y=[lower_row,lower_row+0.3])


    set_title(ax_sharkfin,4,'',(-0.15,0),10)
    axcb1 = plt.axes([0.9,0.1,0.02,0.3])
    axcb1.tick_params(axis='both', which='major', labelsize=12)
    pKL = ax_sharkfin_KL.pcolormesh(res_sharks[options['order'][0]],res_sharks[options['order'][1]],res_sharks['KL_entropy'][p,...],cmap=mpl.cm.get_cmap('Reds',10),shading='auto',vmin=10**(-4),vmax=0.1,norm=mcolors.LogNorm())
    plt.setp(ax_sharkfin_KL,yticklabels=[],xlim=[0,20],xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]')
    set_title(ax_sharkfin_KL,5,'',(-0.15,0),10)

    plt.colorbar(pKL, cax = axcb1)#,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,7))
    axcb1.set_ylabel('KL')

    if save:
        sv_name = './figures/approximations.%s' % (file_format)
        plt.savefig(sv_name,dpi=300)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)

    return res, res_sharks


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
        dI = abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q)))
        min_idx = np.nanargmin(dI)
        min_val = np.nanmin(dI)
        # print(i, nu, min_idx, min_val)
        res['I'][i,0] = np.sqrt(net.I_squared_nu(nu,q[min_idx]))
        res['I'][i,1] = np.sqrt(net.I_squared_nu(nu,nu**2))
        res['q'][i,0] = q[min_idx]
        res['q'][i,1] = net.get_q(nu,nu**2,res['I'][i,0])
        res['alpha'][i,0] = net.alpha(res['q'][i,0])
        res['alpha'][i,1] = net.alpha(res['q'][i,1])

    # fig,axes = plt.subplots(2,2)
    # axes[0,0].plot(res['q'][:,0],'k-')
    # axes[0,0].plot(res['q'][:,1],'k--')
    #
    # axes[1,0].plot(res['alpha'][:,0],'k')
    # axes[1,0].plot(res['alpha'][:,1],'k--')
    # plt.setp(axes[1,0],ylim=[0,0.15])
    #
    # axes[0,1].plot(res['I'][:,0],'k')
    # axes[0,1].plot(res['I'][:,1],'k--')
    # plt.show(block=False)

    return res

def plot_approx(ax_ex,ax_I,ax_alpha,res):
    # print()

    res_approx = get_approximation(0.00)
    res_approx_hetero = get_approximation(0.04)
    # print(res['q_approx'][0,0,:])
    # print(res_approx['q'])
    #
    # print(res['q_approx'][0,2,:])
    # print(res_approx_hetero['q'])

    net = network(0.005,0.01,0.00)
    ylims = [-0.4,0]
    nu = 5
    q = np.linspace(0,4*nu**2,1001)
    y_is = -np.sqrt(net.I_squared_nu(nu,nu**2)) ## intersection with function
    y_is_transformed = (y_is-ylims[0]) / (ylims[1] - ylims[0])
    min_idx = np.nanargmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))

    ax_ex.axvline(nu**2,ymax=y_is_transformed,color=[0.6,0.6,0.6],lw=2,ls='--')
    ax_ex.axhline(-np.sqrt(net.I_squared_nu(nu,nu**2)),color=[0.6,0.6,0.6],lw=2,ls='--')#,label=r"approx.: $\displaystyle q=\bar{\nu}^2$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_nu(nu,q)),'k-',label=r"solution for $\displaystyle \bar{\nu}$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_q(nu,q)),'k--',label=r"solution for $\displaystyle q$")
    ax_ex.annotate(r'$\displaystyle q = \bar{\nu}^2$',xy=[nu**2,-np.sqrt(net.I_squared_nu(nu,nu**2))],xytext=[nu**2-22,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)
    ax_ex.annotate(r'$\displaystyle (\bar{\nu}^{\star},q^{\star})$',xy=[q[min_idx],-np.sqrt(net.I_squared_nu(nu,q[min_idx]))],xytext=[nu**2+15,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)

    plt.setp(ax_ex, ylim=ylims)
    ax_ex.legend(prop={'size':10},bbox_to_anchor=(1.2,0.0),loc='lower right',frameon=False)
    loc = ticker.MultipleLocator(base=0.1) # this locator puts ticks at regular intervals
    ax_ex.yaxis.set_major_locator(loc)
    plot_I(ax_I,res_approx,'k')
    plot_I(ax_I,res_approx_hetero,'r')

    plt.setp(ax_I,xlabel=r'$\displaystyle \bar{\nu}$ [Hz]',ylabel=r'$\displaystyle \frac{\Delta I}{\Psi_0 - I_0}$')

    ax_alpha.plot(res_approx['nu'],res_approx['alpha'][:,0],'k')
    ax_alpha.plot(res_approx['nu'],res_approx['alpha'][:,1],'k--')
    ax_alpha.plot(res_approx['nu'],res_approx_hetero['alpha'][:,0],'r')
    ax_alpha.plot(res_approx['nu'],res_approx_hetero['alpha'][:,1],'r--')


def plot_I(ax,res,ls):
    # nu_lim = 12
    # steps = res['I'].shape[0]
    ax.plot(res['nu'],(res['I'][:,0]-res['I'][:,1])/res['I'][:,0],ls,label=r"$\frac{\Delta I}{I}$")

    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_exact,'k')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,0],'r')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,1],'r--')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),I_approx,'r')
