import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def two_populations(L=1,S=[1,2],steps=100,plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False):

## stats:
####    0: sharkfins
####    1: rate_dependence stats
######    2: borders phase space
####    2: compare exact vs. approx (single)
####    3: KL-phase-space (costly!)

    steps = steps + 1       # correct for removal of first item

    ## general plot setup
    set_plot_params()

    fig,ax = plt.subplots(2,3,figsize=(7.5,7.5),dpi=300)
    plt_para = {
        'ax_label': [],
        'const_label': []
    }

    J_l = np.ones((L,L))
    np.fill_diagonal(J_l,0)
    options = {


        # count of layers, populations, PSPs
        'L': L,
        'S': S,     # contains number of synapses for each population

        # layer level parameters
        'eps': [1./np.sqrt(2),0.3],
        'eta': 0.9,
        'J_l': J_l,
        'kappa': 1.,

        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'alpha_0': 0.02,
        'tau_M': 0.01,
        'J_0': [-1.,1.,-1.,1.],

        # psp level parameters
        'tau_I': [0.03,0.005,0.2,0.03,0.005,0.2],
        'tau_n': 0.1,
        'tau_norm': 1.,

        # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
        'mode_stats': 0,
        'mode_calc': 0,
        'simulation': {
            # for each iteration parameter, specify (layer,population,psp)-tuple
            # specify -1 if a level of hierarchy is non-applicable
            # specify 'None' if should be applied to all candidates
            'tau_n': [0.,1.],
            'sim_prim': [0,1,0],       # when population parameters are iterated, specify population number(s) (empty = all)
            'tau_I': [0.,0.1],
            'sim_sec': [0,1,0],     # when synaptic timeconstants are iterated, specify number within population
        }
    }

    order = list(options['simulation'].keys())

    print(order)


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

    options['rateWnt'] = [2.]
    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    for p in range(2):
        plot_fins(ax[p,1],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)

    options['rateWnt'] = [5.]
    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    for p in range(2):
        plot_fins(ax[p,2],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)


    plt.show(block=False)

    return res
