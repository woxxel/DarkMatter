import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def tau_stats(steps=1000,rateWnt=[2],alpha_0=[0,0.02,0.04],tau_G=[0,0.1],eps=[0.5],eta=[0.9],n=[0],J=-1.,Npop=1,drive=0,save=0,file_format='png',rerun=False,compile=False):

    options = {
        'Npop': Npop,
        'order': ['tau_G','alpha_0','rateWnt','n','eta','eps'],
        'rateWnt': rateWnt,
        'alpha_0': alpha_0,
        'tau_G': tau_G,
        'eps': eps,
        'eta': eta,
        'n': n,
        'drive': drive,
        'mode_stats': 1,
        'J': J
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # return res

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
    print(x_key)

    currents = True

    set_plot_params()
    #mpl.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?

    plt_para = {
        'title': {
            'descr': False,
            'x_offset': -0.2
            },
        'two_row': False
    }

    if plt_para['two_row']:
        fig, ax = plt.subplots(2,3,figsize=(7.5,4),dpi=300)
    else:
        fig, ax = plt.subplots(1,3,figsize=(7.5,1.5),dpi=300)
        # fig.delaxes(ax[1,0])
        # fig.delaxes(ax[1,1])
        # fig.delaxes(ax[1,2])
        # ax[0,0].set_position([0.1,0.25,0.22,0.6])
        # ax[0,1].set_position([0.415,0.25,0.22,0.6])
        # ax[0,2].set_position([0.76,0.25,0.22,0.6])



    #big_ax = fig.add_subplot(111)
    big_ax = plt.axes([0.1,0.13,0.8,0.8])
    big_ax.set_facecolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])

    if x_key == 'rateWnt':
        big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    if x_key == 'tau_G':
        big_ax.set_xlabel(r'$\displaystyle \tau_G\,$[ms]')
    elif x_key == 'eps':
        big_ax.set_xlabel(r'$\displaystyle \varepsilon$')
    elif x_key == 'n':
        big_ax.set_xlabel(r'$\displaystyle b$')

    #print res['trans_implausible']
    #print res
    steps1 = res['gamma'].shape[1]
    pCol = ['r','k']

    ## define dictionary with transition point indices
    trans_idx = {}
    for key in ['inc','imp','DM','np']:
        trans_idx[key] = np.zeros(steps1,'int')
        for a in range(steps1):
            idx = np.where(res[x_key]==res[key+'_trans'][0,a])[0]
            trans_idx[key][a] = idx[0] if len(idx) else -1

    x_lim = options[x_key][-1]
    plot_q(ax[0],res,x_key,trans_idx,plt_para,x_lim,order=0)
    # ax[0,0].text(8,43,r'$\displaystyle q\approx\bar{\nu}^2$',fontsize=10)
    ax[0].legend(prop={'size':10},bbox_to_anchor=(0.4,0.5),loc='lower left',handlelength=1)

    # plot_q_zoom(ax[0,1],res,x_key,trans_idx,plt_para,1,order=0)
    # ax[0,1].text(0.1,3.5,r'$\displaystyle q\approx \frac{\bar{\nu}\nu_{max}}{\sqrt{2}}$',fontsize=10)
    plot_currents(ax[1],res,x_key,trans_idx,plt_para,x_lim,options=['var'],order=0)
    plot_currents(ax[2],res,x_key,trans_idx,plt_para,x_lim,options=['I'],order=0)
    # plot_gamma(ax[2],res,x_key,trans_idx,plt_para,x_lim,order=0)
    # plot_chi(ax[1,1],res,x_key,trans_idx,plt_para,x_lim,order=0)

    set_title(ax[0],1,'',(-0.075,0),10)
    set_title(ax[1],2,'',(-0.075,0),10)
    set_title(ax[2],3,'',(-0.075,0),10)
    # set_title(ax[1,0],4,'',(-0.075,0),10)
    # set_title(ax[1,1],5,'',(-0.075,0),10)
    # set_title(ax[1,2],6,'',(-0.075,0),10)

    plt.setp(ax[0],xlabel='',xticks=np.linspace(0,0.1,3),xticklabels=['%d'%x for x in np.linspace(0,100,3)])
    plt.setp(ax[1],xlabel='',xticks=np.linspace(0,0.1,3),xticklabels=['%d'%x for x in np.linspace(0,100,3)])
    plt.setp(ax[2],xlabel='',xticks=np.linspace(0,0.1,3),xticklabels=['%d'%x for x in np.linspace(0,100,3)])
    # plt.setp(ax[1,0],xlabel='')
    # plt.setp(ax[1,1],xlabel='')

    # options = {
    #     'rateWnt': [0,20],
    #     'alpha_0': [0,0.2],
    #     'tau_G': [0.005],
    #     'J': J
    # }
    # results_bounds = darkMatter(steps=500,options=options,rerun=rerun)
    # plot_regions(ax[1,2],results_bounds,trans_idx,plt_para,res[x_key][-1],order=0)
    # plt.setp(ax[1,2],xlabel='')
    #
    plt.subplots_adjust(left=0.075, bottom=0.25, right=0.95, top=0.8, wspace=0.4, hspace=0.6)
    for j in range(3):
        ax[j].spines['right'].set_color('none')
        ax[j].yaxis.set_ticks_position('left')
        ax[j].spines['top'].set_color('none')
        ax[j].xaxis.set_ticks_position('bottom')

    # fig.tight_layout(pad=0.5)
    if save:
        sv_name = './figures/timescale.%s' % (file_format)
        plt.savefig(sv_name,format=file_format,dpi=300)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)

    return res
