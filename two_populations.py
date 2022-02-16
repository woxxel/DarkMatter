import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def two_populations(steps=100,plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False):

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

    options = {
        'mode_stats': 0,
        'order': ['eps','alpha_0','rateWnt','tau_G','n','eta'],
        'eps': [0,1./np.sqrt(2)],
        'alpha_0': [0,0.1],
        'rateWnt': [1.],
        'tau_G': [0.03],
        'tau_A': [0.005],
        'Npop': 2,
        'eta': [0.9],
        'n': [0],
    }


    results = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    for p in range(2):
        plot_fins(ax[p,0],results[options['order'][0]],results[options['order'][1]],results['gamma'][p,...],results['chi'][p,...],results['regions'][p,...],plt_para)

    # options['rateWnt'] = [2.]
    # results = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # for p in range(2):
    #     plot_fins(ax[p,1],results[options['order'][0]],results[options['order'][1]],results['gamma'][p,...],results['chi'][p,...],results['regions'][p,...],plt_para)
    #
    # options['rateWnt'] = [5.]
    # results = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # for p in range(2):
    #     plot_fins(ax[p,2],results[options['order'][0]],results[options['order'][1]],results['gamma'][p,...],results['chi'][p,...],results['regions'][p,...],plt_para)


    plt.show(block=False)

    return results
