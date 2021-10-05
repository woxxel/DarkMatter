import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.stats as stat
#from scipy.stats import gamma
import matplotlib as mpl


def poisson_distr(x,para):
    
    try:
      return para['lambda']**x*np.exp(-para['lambda'])/math.factorial(x)
    except:
      print 'Please provide integer input for x'
      return


def chi2_distr(x,para):
  #print "under construction"
  
  if x <= 0:
    return 0
  else:
    return (x**(para['n']/2.-1)*np.exp(-x/2.))/(2**(para['n']/2.)*math.gamma(para['n']/2.))
  

def gamma_distr(x,para):
  
  return stat.gamma.pdf(x,para['a'])

def plot_distribution(distr='poisson',para={'lambda':1},plt_range=np.arange(10)):
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    
    plt.figure()
    
    if distr == 'poisson':
      #plt_array = np.linspace(0,10,11)
      func_array = np.zeros(len(plt_range))
      for x in plt_range:
        
        func_array[x] = poisson_distr(x,para)
      
      
    if distr == 'chi2':
        func_array = np.zeros(len(plt_range))
        for i in range(len(plt_range)):
            func_array[i] = chi2_distr(plt_range[i],para)
        
        plt.ylabel = '$\displaystyle \chi^2$'
      
    
    if distr == 'gamma':
      func_array = gamma_distr(plt_range,para)
      
      
    plt.plot(plt_range,func_array)
    plt.show(block=False)