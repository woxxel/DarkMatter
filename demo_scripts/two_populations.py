import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def two_populations(L=1,S=[1,2,2],steps=100,plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False):

## stats:
####    0: sharkfins
####    1: rate_dependence stats
######    2: borders phase space
####    2: compare exact vs. approx (single)
####    3: KL-phase-space (costly!)

    steps = steps       # correct for removal of first item

    ## general plot setup
    set_plot_params()

    fig,ax = plt.subplots(2,3,figsize=(7.5,5),dpi=300)
    plt_para = {
        'ax_label': [],
        'const_label': []
    }

    J_l = np.ones((L,L))
    np.fill_diagonal(J_l,0)

    options = {
        # count of layers, populations, PSPs
        'L': L,
        'P': len(S),
        'S': S,     # contains number of synapses for each population

        # layer level parameters
        'eps': 1./np.sqrt(2),
        'eta': 0.9,
        'J0_l': J_l,

        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'kappa': [1.,0.5,0.5],
        'alpha_0': 0.02,
        # 'Psi_0': [0,0.05,-0.05],
        'Psi_0': [0,0.05],
        'tau_M': 0.01,
        # 'J0': [-1.,1.,1.],
        'J0': [-1.,1.],

        # psp level parameters
        # 'tau_I': [0.01,0.005,0.2,0.005,0.2],
        'tau_I': [0.01,0.005,0.2],
        'tau_n': 0.,
        'tau_norm': 1.,

        # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
        'mode': 0,
        'mode_stats': 0,
        'mode_calc': 0,
        'simulation': {
            # for each iteration parameter, specify (layer,population,psp)-tuple
            # specify -1 if a level of hierarchy is non-applicable
            # specify 'None' if should be applied to all candidates
            'rateWnt': [0.,20.],
            'alpha_0': [0.,0.2],

            'sim_prim': [0,-1,0],       # when population parameters are iterated, specify population number(s) (empty = all)
            'sim_sec': [0,-1,0],     # when synaptic timeconstants are iterated, specify number within population
        }
    }
    print(options)

    order = list(options['simulation'].keys())

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # steps1 = res['gamma'].shape[1]
    #
    # x_key = options['order'][0]
    #
    # fig, ax = plt.subplots(2,3,figsize=(7.5,4),dpi=300)
    #
    # ## define dictionary with transition point indices
    # trans_idx = {}
    # for key in ['inc','imp','DM','np']:
    #     trans_idx[key] = np.zeros(steps1,'int')
    #     for a in range(steps1):
    #         idx = np.where(res[x_key]==res[key+'_trans'][0,a])[0]
    #         trans_idx[key][a] = idx[0] if len(idx) else -1
    #
    # x_lim = res['inc_trans'][0,0]
    # plot_q(ax[0,0],res,x_key,trans_idx,plt_para,x_lim,order=0)
    #
    # plot_currents(ax[0,1],res,x_key,trans_idx,plt_para,x_lim,order=0)
    #
    # plot_gamma(ax[1,0],res,x_key,trans_idx,plt_para,x_lim,order=0)
    # plot_chi(ax[1,1],res,x_key,trans_idx,plt_para,x_lim,order=0)
    #
    # plt.setp(ax[0,0],xlim=options[x_key],ylim=[0,10])
    # plt.setp(ax[0,1],xlim=options[x_key],ylim=[-0.2,0.1])
    # plt.setp(ax[1,0],xlim=options[x_key])
    # plt.setp(ax[1,1],xlim=options[x_key],ylim=[0,10])

    # for i in range(res['gamma'].shape)
    for p in range(2):
        plot_fins(ax[p,0],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)

    # options['tau_I'] = [0.03,0.005,0.2]
    # res = darkMatter(steps=steps,options=options,rerun=rerun,compile=False)
    # for p in range(2):
    #     plot_fins(ax[p,1],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)
    #
    # options['tau_I'] = [0.06,0.005,0.2]
    # res = darkMatter(steps=steps,options=options,rerun=rerun,compile=False)
    # for p in range(2):
    #     plot_fins(ax[p,2],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)


    big_ax = fig.add_axes([0.1,0.1,0.8,0.85])
    big_ax.set_facecolor('none')
    big_ax.tick_params(labelcolor='none',top=False,bottom=False,left=False,right=False)
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)

    # plt.setp(big_ax,xlabel=r'$\displaystyle \bar{\nu}$',ylabel=r'$\displaystyle r$')
    plt.setp(big_ax,xlabel=r"$\bar{ \varepsilon }$",ylabel=r'$\alpha_0$')

    fig.tight_layout(pad=3.0)
    if save:
        sv_name = './figures/two_pop_n.%s' % (file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)

    return res
