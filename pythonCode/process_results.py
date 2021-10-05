import numpy as np
#import sympy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import scipy.stats as stats
import imp, math, time

from operator import itemgetter, attrgetter, methodcaller

imp.load_source('read_mat', 'read_mat.py')
from read_mat import *

def get_data(fileNames):
    
    # get data into proper shape
    
    N = len(fileNames)
    print N
    
    data = {}
    data['dt'] = np.zeros(N)
    #data['n_neuron'] = np.zeros(N).astype('int')
    data['n_frame'] = np.zeros(N).astype('int')
    
    data['ROI_label'] = [[] for i in range(N)]
    data['events'] = [[] for i in range(N)]

    idx = 0
    for dat in fileNames:
        print dat
        
        data_tmp = read_mat(dat)
        data['dt'][idx] = 1./data_tmp['frame_rate'][0]
        
        n_neuron, data['n_frame'][idx] = data_tmp['data'][0,0].shape
        
        data['events'][idx] = [[] for i in range(n_neuron)]
        data['ROI_label'][idx] = list(data_tmp['ROI_label'][:,0]-1)
        
        for n in range(n_neuron):
            
            event_tmp = data_tmp['events'][0,n].reshape((2,-1))
            
            data['events'][idx][n].extend(event_tmp[0]*data['dt'][idx])
            
        idx += 1
    
    return data


def stretch(x,k):   # make moebius trafo (ax+b)/(cx+d) to transform to intervall [0,1), therefore T(0) = 0 -> b=0, T(infty) = a/c -> a=c, a/d=k
    
    
    return x/(x+k)


def get_lambda(n,T,alpha=None):
  
  mean_lambda = n/float(T)
  
  if alpha:
    if not n:
      left_lambda = 0
    else:
      left_lambda = 1/2.*stats.chi2.interval(1-alpha,2*n)[0]/float(T)
    
    right_lambda = 1/2.*stats.chi2.interval(1-alpha,2*(n+1))[1]/float(T)
  
  else:
    left_lambda = 0
    right_lambda = 0
  
  return [mean_lambda,left_lambda,right_lambda]


def plt_CI(X,Y,ax,col='k',ls='-',ms='o',label=None,mask=None):
  
  assert len(X.shape) == 1, 'Please specify a one dimensional range array!'
  assert Y.shape[0] == X.shape[0], 'X and Y arrays should have the same length!'
  assert Y.shape[1] == 3, 'Y array should include mean value and upper and lower bound of confidence intervall!'

  if mask == None:
    mask = np.ones(len(X)).astype('bool')

  ax.plot(X[mask],Y[mask,0],'-',color=col,linestyle=ls,marker=ms,label=label,linewidth=2)
  
  ax.fill_between(X[mask],Y[mask,1],Y[mask,2],color=col,linestyle=ls,alpha=0.2,edgecolor=None,linewidth=0)
  
  #if fit_func:
    #Y_sigma = np.max(Y[:,1:] - Y[:,0].reshape(len(Y),1),axis=1)
    
    #popt,pcov = curve_fit(fit_func,X[mask],Y[mask,0],sigma=Y_sigma[mask],p0=p0)
    
    #perr = np.sqrt(np.diag(pcov))
    
    #print 'fit results: ', popt
    #print 'fit errors: ', perr
    #ax.plot(X,fit_func(X,popt[0],popt[1]),'--',color='r',linewidth=0.5)
    
    #return popt,perr


def plot_frame(ax,hide={'spines':1,'ticks':1}):
    """Hides the top and rightmost axis spines from view for all active
    figures and their respective axes."""

    # Retrieve a list of all current figures.
    #figures = [x for x in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    #for figure in figures:
      # Get all Axis instances related to the figure.
      #for ax in figure.canvas.figure.get_axes():
        # Disable spines.
    if hide['spines']:
      ax.spines['right'].set_color('none')
      ax.spines['top'].set_color('none')
    if hide['spines'] == 2:
      ax.spines['left'].set_color('none')
      ax.spines['bottom'].set_color('none')
    
    if hide['ticks']:
    # Disable ticks.
      ax.xaxis.set_ticks_position('bottom')
      ax.yaxis.set_ticks_position('left')
    if hide['ticks'] == 2:
      ax.xaxis.set_ticks_position('none')
      ax.yaxis.set_ticks_position('none')
      
            
  
def process_results(fileNames,plot=1,save=0):
    
    N = len(fileNames)
    
    data = get_data(fileNames)
    
    #return data
    data['n_neuron'] = 0
    data['t_offset'] = np.zeros(N)

    print data['dt']
    print data['n_frame']
    
    for idx in range(N):
      
      data['n_neuron'] = max(data['n_neuron'],data['ROI_label'][idx][-1] + 1)
      data['t_offset'][idx] = data['dt'][idx]*data['n_frame'][idx]
    
      #data['t_offset'][idx] = 0
      #for i in range(idx+1):
        
    print 'offset: ', data['t_offset']
    print 'neurons: ', data['n_neuron']
    
    
    data['spikeTrain'] = [[] for i in range(data['n_neuron'])]
    data['N_ap'] = np.zeros((N,data['n_neuron']))
    data['firing_rate'] = np.zeros((N,data['n_neuron'],3))
    #data['N_ap_init'] = np.zeros(data['n_neuron'])
    t_border = 10000 #data['t_offset'][1]
    
    t_max = 0
    for i in range(N):
      t_max += data['dt'][i]*data['n_frame'][i]
    print "t_max: ", t_max
    
    #print N
    for idx in range(N):
      for n in range(len(data['ROI_label'][idx])):
        data['spikeTrain'][data['ROI_label'][idx][n]].extend([x+idx*data['t_offset'][idx] for x in data['events'][idx][n]])
        #print data['spikeTrain'][data['ROI_label'][idx][n]]
        #print sum(data['spikeTrain'][data['ROI_label'][idx][n]] < t_border)
        #if idx == 0:
          #data['N_ap_init'][data['ROI_label'][idx][n]] = len(data['spikeTrain'][data['ROI_label'][idx][n]])
        #elif idx == N-1:
        data['N_ap'][idx,data['ROI_label'][idx][n]] = len(data['events'][idx][n])
        #data['N_ap'][idx,data['ROI_label'][idx][n]] = len(data['spikeTrain'][data['ROI_label'][idx][n]])
        #if idx==1:
            #print data['N_ap'][idx,data['ROI_label'][idx][n]]
            #data['N_ap'][idx,data['ROI_label'][idx][n]] -= data['N_ap'][0,data['ROI_label'][idx][n]]
    
    #print data['spikeTrain']
    # now get observables: spike count correlation(?!), histogram of ISI, ...
    
    # spike count correlation
    num_bins = int(t_max)*1
    bins = np.linspace(0,int(t_max),num_bins)
    data['bin_count'] = np.zeros((data['n_neuron'],num_bins))
    
    data['corr'] = np.zeros((data['n_neuron'],data['n_neuron']))
    for n in range(data['n_neuron']):
      for e in data['spikeTrain'][n]: 
        data['bin_count'][n,np.where(e > bins)[0][-1]] += 1
    
    for n_ref in range(data['n_neuron']):
      for n in range(n_ref,data['n_neuron']):
        data['corr'][n_ref,n] = np.corrcoef(data['bin_count'][n_ref],data['bin_count'][n])[0,1]
        if np.isnan(data['corr'][n_ref,n]):
          data['corr'][n_ref,n] = 0
        data['corr'][n,n_ref] = data['corr'][n_ref,n]
    

    # sort stuff by correlation population
    data['corr_population'] = [[]]*data['n_neuron']
    for n in range(data['n_neuron']):
      data['corr_population'][n] = (n,np.sum(data['corr'][n] > 0.5),np.nansum(data['corr'][n]))
    
    list_tmp = sorted(data['corr_population'],key=itemgetter(1,2))
    
    new_index = np.zeros((data['n_neuron'],2)).astype('int')
    for n in range(data['n_neuron']):
      new_index[n] = (list_tmp[n][0],n)
    
    data['corr_sort'] = np.zeros((data['n_neuron'],data['n_neuron']))
    for n_ref in range(data['n_neuron']):
      for n in range(data['n_neuron']):
        #print new_index[n_ref]
        data['corr_sort'][n_ref,n] = data['corr'][new_index[n_ref,0],new_index[n,0]]
    
    
    #print new_index
    
    min_val = -0.2
    max_val = 0.7
    
    zero_val = (0-min_val)/(max_val-min_val)
      
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    
    
    if plot == 1:
        levs = range(160)
        rb = mcolors.LinearSegmentedColormap.from_list(name='red_blue',colors=[(0,(0,0,1)),(1,(1,0,0))],N=len(levs)-1,)
        
        rwb = mcolors.LinearSegmentedColormap.from_list(name='red_white_blue',colors=[(0,(0,0,1)),(zero_val,(1,1,1)),(1,(1,0,0))],N=len(levs)-1,)
        
        plt.figure(figsize=(8,6))
        
        ax1 = plt.axes([0.1,0.1,0.35,0.85])
        ax2 = plt.axes([0.65,0.55,0.3,0.4])
        ax22 = plt.axes([0.81,0.775,0.14,0.15])
        ax3 = plt.axes([0.65,0.1,0.3,0.35])
        #ax4 = plt.axes([0.775,0.1,0.2,0.4])
        #ax5 = plt.axes([0.775,0.55,0.2,0.4])
        
        ax_cb = plt.axes([0.475,0.1,0.025,0.85])
        #ax_cb2 = plt.axes([0.7,0.1,0.25,0.025])
        
        if N > 1:
            offset = data['t_offset'][1]
        else:
            offset = t_max
        
        alpha = 0.05
        
        for n in range(data['n_neuron']):
        
            for i in range(N):
                    data['firing_rate'][i,n] = get_lambda(data['N_ap'][i,new_index[n,0]],data['t_offset'][i],0.05)
            
            scale = 3.
            c = (max(np.log10(data['firing_rate'][0,n,0]),-scale) + scale)/scale
            #print c
            col = (c,0,1-c)
            
            ax1.plot(data['spikeTrain'][new_index[n,0]],n*np.ones(np.sum(data['N_ap'][:,new_index[n,0]],axis=0)),'|',color=col,markersize=5)
        
        ax1.plot([offset, offset],[0,data['n_neuron']],'k--')
        ax1.set_xlim([0, t_max])
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('\# neuron')
        plot_frame(ax1,hide={'ticks':2,'spines':0})
        
        norm = mpl.colors.Normalize(vmin=-scale, vmax=0)
        cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=rb,norm=norm,orientation='vertical')
        cb1.set_label(r'$\displaystyle \log_{10}(\nu_1)$')
        
        
        test = np.array(sorted(data['firing_rate'][0],key=itemgetter(0)))
        #print test
        
        #plt.figure()
        #ax22 = plt.axes([0.125,0.15,0.7,0.8])
        plt_CI(np.arange(data['n_neuron']),test,ax22,ms=None,label='firing rate')
        ax22.set_xlabel(r'\# neuron')
        ax22.set_ylabel(r'$\displaystyle \nu_1$ ($\displaystyle s^{-1}$)')
        ax22.set_xticks(np.linspace(0,80,5))
        ax22.set_yticks(np.linspace(0,0.15,4))
        #plt.show(block=False)
        
        
        #plt_CI(np.arange(data['n_neuron']),data['firing_rate'][0],ax22,ms=None,label='firing rate')
        #ax22.set_xlabel(r'\# neuron')
        #ax22.set_ylabel(r'$\displaystyle \nu_1$ ($\displaystyle s^{-1}$)')
        #ax22.set_xlim([0,data['n_neuron']])
        #ax22.set_xticks(np.linspace(0,80,5))
        #ax22.set_yticks(np.linspace(0,0.15,4))
        
        plot_frame(ax22,hide={'ticks':1,'spines':1})

        #mask = np.isfinite(np.log(data['firing_rate']))
        ax2.hist(data['firing_rate'][0,:,0],bins=np.linspace(0,0.2,101),color='r',alpha=0.6,label=r'$\displaystyle \nu_1$')
        if N>1:
            ax2.hist(data['firing_rate'][1,:,0],bins=np.linspace(0,0.2,101),color='b',alpha=0.6,label=r'$\displaystyle \nu_2$')
        ax2.set_xlim([0,0.15])
        ax2.set_xticks(np.linspace(0,0.15,4))
        ax2.set_yticks(np.linspace(0,data['n_neuron']/2.,6))
        ax2.set_yticklabels(np.linspace(0,0.5,6))
        ax2.set_xlabel(r'$\displaystyle \nu$ (Hz)')
        ax2.set_ylabel(r'$\displaystyle \rho(\nu)$')
        ax2.legend(loc=4,prop={'size':12})
        plot_frame(ax2)
        
        
        if N>1:
            
            scale = 3.
            
            for n in range(data['n_neuron']):
            
                    c = (max(np.log10(data['firing_rate'][0,n,0]),-scale) + scale)/scale
                    col = (c,0,1-c)
                    
                    ax3.plot(data['firing_rate'][0,n,0],data['firing_rate'][1,n,0],'o',color=col,markersize=5)
            
            x_err = np.abs([data['firing_rate'][0,:,1] - data['firing_rate'][0,:,0],data['firing_rate'][0,:,2] - data['firing_rate'][0,:,0]])
            y_err = np.abs([data['firing_rate'][1,:,1] - data['firing_rate'][1,:,0],data['firing_rate'][1,:,2] - data['firing_rate'][1,:,0]])
            
            ax3.errorbar(data['firing_rate'][0,:,0],data['firing_rate'][1,:,0],xerr=x_err,yerr=y_err,ecolor='k',fmt='none')
            
            ax3.set_xlabel(r'$\displaystyle \nu_1$ (Hz)')
            ax3.set_ylabel(r'$\displaystyle \nu_2$ (Hz)')
            #ax3.text(0.06,0.04,r'$\displaystyle n(\nu_1=0)$: %d'%np.sum(data['firing_rate'][0,:,0]==0))
            #ax3.text(0.06,0.02,r'$\displaystyle n(\nu_2=0)$: %d'%np.sum(data['firing_rate'][1,:,0]==0))
            
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            
            ax3.set_xlim([10**(-3),0.3])
            ax3.set_ylim([10**(-3),0.3])
        
        #cbar0 = plt.colorbar(cbar_dummy)
        #plt.colorbar(ax=ax1)
        if save==1:
            sv_name = './../pics/spike_analysis.pdf'
            plt.savefig(sv_name)
            print "Figure saved as %s" % sv_name
        
        plt.show(block=False)
        
        plt.figure(figsize=(4,3))
        ax = plt.axes([0.125,0.15,0.7,0.8])
        ax_cb = plt.axes([0.8,0.15,0.025,0.8])
        ax.imshow(data['corr'],cmap=rwb,vmin=-1,vmax=1)
        ax.imshow(data['corr_sort'],cmap=rwb,vmin=min_val,vmax=max_val,interpolation='none')
        ax.set_xlabel(r'\# neuron')
        ax.set_ylabel(r'\# neuron')
        norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
        cb = mpl.colorbar.ColorbarBase(ax_cb, cmap=rwb,norm=norm,orientation='vertical')
        #cb2.set_label(r'$\displaystyle \log_{10}(\nu)$')
        cb.set_label('Pearson Correlation')
        cb.set_ticks(np.linspace(-1,1,11))
        
        if save==2:
            sv_name = './../pics/correlation.pdf'
            plt.savefig(sv_name)
            print "Figure saved as %s" % sv_name
        
        plt.show(block=False)
    
    
    ## construct cummulative density function
    #print "N_ap: ", 
    
    N_ap = np.array(sorted(np.sum(data['N_ap'],axis=0)))
    
    #print "rates: ", 
    
    rates = sorted(np.sum(data['N_ap'],axis=0))/np.sum(data['t_offset'])
    max_N_ap = int(max(np.sum(data['N_ap'],axis=0)))
    
    N_ap_hist = np.histogram(np.sum(data['N_ap'],axis=0),bins=np.arange(max_N_ap+2))
    
    data['firing_rate'] = np.zeros((max_N_ap+1,4))
    
    for n in range(max_N_ap+1):
        data['firing_rate'][n,0] = N_ap_hist[0][n]
        data['firing_rate'][n,1:] = get_lambda(N_ap_hist[1][n],np.sum(data['t_offset']),0.05)
    
    x_max = 1.      # value, at which cdf = 1. How does it influence the solution?
    #cdf = np.zeros(data['n_neuron'])
    p_array = np.cumsum(data['firing_rate'][:,0]/float(data['n_neuron']+1))
    
    print data['firing_rate']
    data['cdf'] = np.zeros((max_N_ap+2,3))
    data['cdf'][:-1,0] = p_array                                                                        # p
    data['cdf'][:-1,1] = data['firing_rate'][:,1]                                                       # rate
    data['cdf'][:-1,2] = 1./(data['n_neuron']+1)*np.sqrt((data['n_neuron']*p_array*(1-p_array)))        # sigma
    #print p_array*(1-p_array)
    # set boundary conditions
    data['cdf'][-1,0] = 1
    data['cdf'][-1,1] = x_max
    data['cdf'][-1,2] = 0
    
    print p_array
    print data['cdf']
    #print [
    cdf_cut = data['cdf'][N_ap.astype('int')]
    
    #print cdf_cut.shape
    
    plt.figure()
    plt.plot(data['cdf'][:,1],data['cdf'][:,0])
    plt.show(block=False)
    #return data, cdf_cut
    
    #return data
    # define matrices A, B and vector c
    
    #A = sp.zeros(data['n_neuron'],data['n_neuron'])p_array
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
    return A, B, c#, inv
    #return data
    
    #test = np.array(sorted(data['firing_rate'][0],key=itemgetter(0)))
    
    #plt.figure()
    #ax = plt.axes([0.125,0.15,0.7,0.8])
    #plt_CI(np.arange(data['n_neuron']),test,ax,ms=None,label='firing rate')
    #ax.set_xlabel(r'\# neuron')
    #ax.set_ylabel(r'$\displaystyle \nu_1$ ($\displaystyle s^{-1}$)')
    #ax.set_xticks(np.linspace(0,80,5))
    #ax.set_yticks(np.linspace(0,0.15,4))
    #plt.show(block=False)
    #
    
    #plt.figure(figsize=(6,4))
    #ax = plt.axes([0.15,0.15,0.8,0.8])
    #if N>1:
        
      #scale = 3.
      
      #for n in range(data['n_neuron']):
    
        #c = (max(np.log10(data['firing_rate'][0,n,0]),-scale) + scale)/scale
        
        #col = (c,0,1-c)
        
        #ax.plot(data['firing_rate'][0,n,0],data['firing_rate'][1,n,0],'o',color=col,markersize=5)
      
      #x_err = np.abs([data['firing_rate'][0,:,1] - data['firing_rate'][0,:,0],data['firing_rate'][0,:,2] - data['firing_rate'][0,:,0]])
      #y_err = np.abs([data['firing_rate'][1,:,1] - data['firing_rate'][1,:,0],data['firing_rate'][1,:,2] - data['firing_rate'][1,:,0]])
      
      #ax.errorbar(data['firing_rate'][0,:,0],data['firing_rate'][1,:,0],xerr=x_err,yerr=y_err,ecolor='k',fmt='none')
      
      #ax.set_xlabel(r'$\displaystyle \nu_1$')
      #ax.set_ylabel(r'$\displaystyle \nu_2$')
      #ax.text(0.06,0.04,r'$\displaystyle n(\nu_1=0)$: %d'%np.sum(data['firing_rate'][0,:,0]==0))
      #ax.text(0.06,0.02,r'$\displaystyle n(\nu_2=0)$: %d'%np.sum(data['firing_rate'][1,:,0]==0))
      
      #ax.set_xscale('log')
      #ax.set_yscale('log')
      
      #ax.set_xlim([10**(-3),1])
      #ax.set_ylim([10**(-3),1])
      
      #plt.show(block=False)
      
      
    #return data
        #for n in range(data['n_neuron']):
        #spikeTrain = 