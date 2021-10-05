import numpy as np
import scipy as sp
import scipy.misc
import math
import matplotlib.pyplot as plt

def poisson(lam,k):   # only single values allowed
    
    p = np.exp(k*np.log(lam)-lam-math.lgamma(k+1)) # numerically stable computation (avoiding overflow from power and factorial)
    
    return p


def plot_poisson(y):
    
    x_max = 50
    x_array = np.arange(x_max)
    
    p_lam_array = np.zeros(x_max)
    p_k_array = np.zeros(x_max)
    for x in x_array:
        print x
        p_lam_array[x] = poisson(y,x)
        p_k_array[x] = poisson(x,y)
    
    plt.figure()
    plt.plot(x_array,p_lam_array,color='k',label='p_lam')
    plt.plot(x_array,p_k_array,color='r',label='p_k')
    plt.legend()
    plt.show(block=False)