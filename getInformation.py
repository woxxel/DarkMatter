import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *

def information(steps=100,rateWnt=[0,20],alpha_0=[0],tau_G=[0.005],eps=[0.5],eta=[0.9],n=[0],J=-1.,Npop=1,drive=0,save=0,file_format='png',rerun=False,compile=False):

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
        'mode_stats': 4,
        'J': J
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    return res
