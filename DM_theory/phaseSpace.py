import numpy as np
from .functions import *

def create_phaseSpace(steps=201,gamma_range=[0,2],delta_range=[0,6],nu_max=20.):
    ## delta / gamma phase space
    delta_x, gamma_y = np.meshgrid(np.linspace(*delta_range,steps), np.linspace(*gamma_range,steps), indexing='xy')
    chi = get_chi(gamma_y,delta_x,nu_max)
    return chi
