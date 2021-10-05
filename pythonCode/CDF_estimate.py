import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import scipy.stats as stats
import imp, math, time

def get_estimate(alpha_0=0.0,N=100,T=100,n_bin=201):
    
    data = get_samples_from_theory(alpha_0=alpha_0,N=N,T=T,n_bin=n_bin)
    
    CDF = np.cumsum(np.ones(len(data['rate_inf'])))/float(len(data['rate_inf']))
    
    plt.figure()
    plt.plot(sorted(data['rate_inf']),CDF)
    plt.show(block=False)
    
    return
    
    #N_ap = np.array(sorted(np.sum(data['N_ap'],axis=0)))
    
    #print "rates: ", 
    
    #rates = sorted(np.sum(data['N_ap'],axis=0))/np.sum(data['t_offset'])
    #max_N_ap = int(max(np.sum(data['N_ap'],axis=0)))
    
    #N_ap_hist = np.histogram(np.sum(data['N_ap'],axis=0),bins=np.arange(max_N_ap+2))
    
    #data['firing_rate'] = np.zeros((max_N_ap+1,4))
    
    #for n in range(max_N_ap+1):
        #data['firing_rate'][n,0] = N_ap_hist[0][n]
        #data['firing_rate'][n,1:] = get_lambda(N_ap_hist[1][n],np.sum(data['t_offset']),0.05)
    
    x_max = 1.      # value, at which cdf = 1. How does it influence the solution?
    #cdf = np.zeros(data['n_neuron'])
    p_array = np.cumsum(data['firing_rate'][:,0]/float(data['n_neuron']+1))
    
    data['cdf'] = np.zeros((max_N_ap+2,3))
    data['cdf'][:-1,0] = p_array                                                                        # p
    data['cdf'][:-1,1] = data['firing_rate'][:,1]                                                       # rate
    data['cdf'][:-1,2] = 1./(data['n_neuron']+1)*np.sqrt((data['n_neuron']*p_array*(1-p_array)))        # sigma
    print p_array*(1-p_array)
    # set boundary conditions
    data['cdf'][-1,0] = 1
    data['cdf'][-1,1] = x_max
    data['cdf'][-1,2] = 0
    
    
    #print [
    cdf_cut = data['cdf'][N_ap.astype('int')]
    
    #print cdf_cut.shape
    
    plt.figure()
    plt.plot(data['cdf'][:,1],data['cdf'][:,0])
    plt.show(block=False)
    #return data, cdf_cut
    
    #return data
    # define matrices A, B and vector c
    
    #A = sp.zeros(data['n_neuron'],data['n_neuron'])
    #B = sp.zeros(data['n_neuron'],data['n_neuron'])
    #c = sp.zeros(data['n_neuron'],1)
    
    #z = sp.symbols('lambda')
    
    A = np.zeros((data['n_neuron'],data['n_neuron']))
    B = np.zeros((data['n_neuron'],data['n_neuron']))
    c = np.zeros(data['n_neuron'])
    
    for l in range(1,data['n_neuron']+1):
        for m in range(1,data['n_neuron']+1):
            # construct A
            factor = 2*math.pi**4*l**2*m**2
            #print factor
            if l == m:
                A[l-1,m-1] = factor*(1./2-np.sin(2*m*math.pi)/(4.*m*math.pi))
            else:
                A[l-1,m-1] = factor*(1./(2*math.pi*(l**2-m**2)))*((l+m)*np.sin((l-m)*math.pi)-(l-m)*np.sin((l+m)*math.pi))
            
            # construct B
            B[l-1,m-1] = 2./data['n_neuron']*np.sum(np.sin(m*math.pi*cdf_cut[:,1])*np.sin(l*math.pi*cdf_cut[:,1])/cdf_cut[:,2]**2)
            
        # construct c
        tmp = 2./data['n_neuron']*np.sum((cdf_cut[:,0]-cdf_cut[:,1])*np.sin(l*math.pi*cdf_cut[:,1])/cdf_cut[:,2]**2)
        #print tmp
        c[l-1] = 2./data['n_neuron']*np.sum((cdf_cut[:,0]-cdf_cut[:,1])*np.sin(l*math.pi*cdf_cut[:,1])/cdf_cut[:,2]**2)
    
    #inv = z*(A-z*B)**-1*c
    #print inv
    #plt.figure()
    #plt.imshow(A)
    #plt.show()
    return A, B, c