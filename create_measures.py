import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def create_measures(L=1,S=[1,2],plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False):

    steps = 1

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
        'eps': 1./np.sqrt(2),
        'eta': 0.9,
        'J0_l': J_l,
        'kappa': 1.,

        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'alpha_0': 0.02,
        'tau_M': 0.01,
        'J0': 1.,

        # psp level parameters
        'tau_I': [0.01,0.005,0.2],
        'tau_n': 0.,
        'tau_norm': 1.,

        # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
        'mode': 1,      # mode 0: phase plots, mode 1: data simulation
        'mode_stats': 0,
        'mode_calc': 0,
        'simulation': {

            # 'seed_theory': 1,
            # 'seed_time': 1.,
        },
        'computation': {
            'N': 50,    # number of neurons to draw from
            'T': 600.,   # time of measurement

            'draw_from_theory': 4,
            'draw_finite_time': 1,
        }
    }

    # order = list(options['simulation'].keys())
    #
    # print(order)


    res = darkMatter(steps=steps,options=options,mode=1,rerun=rerun,compile=compile)
    return res
