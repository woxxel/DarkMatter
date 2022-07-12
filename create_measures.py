import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def create_measures(L=1,S=[1,2],N=100,T=600.,num_animals=4,plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False,**opts):

    steps = 1

    print(opts)

    ## general plot setup
    # set_plot_params()
    #
    # fig,ax = plt.subplots(2,3,figsize=(7.5,7.5),dpi=300)
    # plt_para = {
    #     'ax_label': [],
    #     'const_label': []
    # }

    J_l = np.ones((L,L))
    np.fill_diagonal(J_l,0)

    options = {
        # count of layers, populations, PSPs
        'L': L,
        'P': 2,
        'S': S,     # contains number of synapses for each population

        # layer level parameters
        'eps': opts['eps'] if 'eps' in opts else 0., #1./np.sqrt(2),
        'eta': opts['eta'] if 'eta' in opts else 0.9,
        'J0_l': J_l,
        'kappa': opts['kappa'] if 'kappa' in opts else 1.,

        # population level parameters
        'I_ext': 1,
        'rateWnt': opts['rateWnt'] if 'rateWnt' in opts else 1.,
        'alpha_0': opts['alpha_0'] if 'alpha_0' in opts else 0.02,
        'tau_M': opts['tau_M'] if 'tau_M' in opts else 0.01,
        'J0': 1.,

        # psp level parameters
        'tau_I': opts['tau_I'] if 'tau_I' in opts else [0.01,0.005,0.2],
        'tau_n': opts['tau_n'] if 'tau_n' in opts else 0.,
        'tau_norm': opts['tau_norm'] if 'tau_norm' in opts else 1.,

        # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
        'mode': 1,      # mode 0: phase plots, mode 1: data simulation
        'mode_stats': 0,
        'mode_calc': 0,
        'simulation': {

            # 'seed_theory': 1,
            # 'seed_time': 1.,
        },
        'computation': {
            'N': N,    # number of neurons to draw from
            'T': T,   # time of measurement

            'draw_from_theory': num_animals,
            'draw_finite_time': 1,
        }
    }
    print(options)

    # order = list(options['simulation'].keys())
    #
    # print(order)


    res = darkMatter(steps=steps,options=options,mode=1,rerun=rerun,compile=compile)
    return res
