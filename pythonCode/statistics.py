import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib as mpl

def plt_CI(X,Y,ax,col='k',ls='-',lw=1,ms='o',mSz=4,label=None,mask=None):
  
  assert len(X.shape) == 1, 'Please specify a one dimensional range array!'
  assert Y.shape[0] == X.shape[0], 'X and Y arrays should have the same length!'
  assert Y.shape[1] == 3, 'Y array should include mean value and upper and lower bound of confidence intervall!'

  if mask == None:
    mask = np.ones(len(X)).astype('bool')

  ax.plot(X[mask],Y[mask,0],'-',color=col,linestyle=ls,linewidth=lw,marker=ms,markersize=mSz,label=label)
  
  ax.fill_between(X[mask],Y[mask,1],Y[mask,2],color=col,linestyle=ls,alpha=0.2,edgecolor=None,linewidth=0)
  
  
  
def poisson_CI(N_max,T=1,mode='T',save=0):
  
    alpha = 0.05            # confidence level
    N_ap = np.arange(N_max)   # number of spikes in interval of length T
    
    T = float(T)
    rate = N_ap/T
    CI = np.array([rate,0.5*stats.chi2.interval(1-alpha,2*N_ap)[0]/T,0.5*stats.chi2.interval(1-alpha,2*(N_ap+1))[1]/T]).transpose() ###### check that
    
    if mode=='var':
        CI_normal = np.array([rate,(N_ap-1.96*np.sqrt(N_ap))/T,(N_ap+1.96*np.sqrt(N_ap))/T]).transpose()
        CI_rel_normal = (CI_normal[:,2]-CI_normal[:,1])/rate
    elif mode=='T':
        rate_ref = N_ap
        CI_ref = np.array([rate_ref,0.5*stats.chi2.interval(1-alpha,2*N_ap)[0],0.5*stats.chi2.interval(1-alpha,2*(N_ap+1))[1]]).transpose()
        CI_rel_ref = (CI_ref[:,2]-CI_ref[:,1])/rate_ref
    elif mode=='T_scale':
        N_ap_ref_sc = N_ap*10
        rate_ref_sc = N_ap_ref_sc/(10.*T)
        CI_ref_sc = np.array([rate_ref_sc,0.5*stats.chi2.interval(1-alpha,2*N_ap_ref_sc)[0]/(10.*T),0.5*stats.chi2.interval(1-alpha,2*(N_ap_ref_sc+1))[1]/(10.*T)]).transpose()
        CI_rel_ref_sc = (CI_ref_sc[:,2]-CI_ref_sc[:,1])/rate_ref_sc
        
    #print N_ap_ref_sc
    #print rate_ref_sc
    CI_rel = (CI[:,2]-CI[:,1])/rate
    
    # plot stuff
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    
    plt.figure(figsize=(3,3))
    ax = plt.axes([0.225,0.575,0.75,0.4])
    ax1 = plt.axes([0.225,0.15,0.75,0.4])
    
    if mode=='var':
        plt_CI(rate,CI_normal,ax,col='r')
    elif mode=='T':
        plt_CI(rate_ref,CI_ref,ax,col='r')
    elif mode=='T_scale':
        plt_CI(rate_ref_sc,CI_ref_sc,ax,col='r')
        
    plt_CI(rate,CI,ax,col='k')
    ax.set_xticks([])
    ax.set_xlim([0,rate[-1]])
    ax.set_ylim([0,CI[-1,0]*1.1])
    #ax.set_xlabel(r'$\displaystyle \nu$')
    ax.set_ylabel(r'$\displaystyle \nu\pm CI$ (Hz)')
    #plt.show(block=False)
    
    if mode=='var':
        ax1.plot(rate,CI_rel_normal,'r',label=r'$\displaystyle T=%g (1.96\sigma)$'%T)
    elif mode=='T':
        ax1.plot(rate_ref,CI_rel_ref,'r',label=r'$\displaystyle T=1$')
    elif mode=='T_scale':
        ax1.plot(rate_ref_sc,CI_rel_ref_sc,'r',label=r'$\displaystyle T=%g,\, (N=%d)$'%(T*10,N_max*10))
    
    ax1.plot(rate,CI_rel,'k',label=r'$\displaystyle T=%g$'%T)
    ax1.set_xlim([0,rate[-1]])
    ax1.set_xlabel(r'$\displaystyle \nu$ (Hz)')
    ax1.set_ylabel(r'$\displaystyle \frac{\Delta CI}{\nu}$')
    ax1.legend(loc=1,prop={'size':12})
    
    if save==1:
        sv_name = './../pics/poisson_CI_%s.pdf'%mode
        plt.savefig(sv_name)
        print "Figure saved as %s" % sv_name
        
    plt.show(block=False)