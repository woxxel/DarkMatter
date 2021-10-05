import math
import matplotlib.pyplot as plt

def DM_cdf(rate=1,steps=3,lower=0.1,upper=0.8):
    #steps = 6
    #rate = np.linspace(0.5,2,steps)
    alpha_0 = np.linspace(0,0.15,steps)#np.linspace(0.02,0.02,steps)#
    tau_G = np.linspace(5,5,steps)#np.linspace(5,60,steps)#
    
    array = alpha_0
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
  
    fig, ax = plt.subplots(2,1,figsize=(8,5))
    fig.delaxes(ax[1])
    axhyper = plt.axes([0.55,0.1,0.35,0.3])
    axsilent = plt.axes([0.1,0.1,0.35,0.3])
    
    ax[0].set_position([0.1,0.5,0.8,0.4])
    silent = np.zeros(steps)
    hyper = np.zeros(steps)
    
    
    
    for i in range(steps):
        p_range, p_exact, cdf = sharkfins(steps=1,rateWnt=[[rate]],alpha_0=[[alpha_0[i]]],tau_G=[[tau_G[i]]],n=[[0]],eps=[[0]],eta=[[0]],tau_M=10,Npop=1,mode_calc='exact',mode_stats=3,save=0,file_format='png')
        
        nu_max = 1./(2*math.pi*math.sqrt(tau_G[i]/1000.*0.010))
        print "nu max: ", nu_max
        
        #print p_range
        d_nu = p_range[1]-p_range[0]
        ### resize stuff
        p_exact = p_exact/cdf[-1]
        cdf = cdf/max(cdf)
        col = i/float(steps)
        
        ax[0].plot(p_range,1-cdf,color=(col,col,col),label=r"$\displaystyle \alpha_0 = %4.2f$"%array[i])
        ax[0].set_yscale('log')
        ax[0].set_ylim([10**(-3),10**0])
        
        idx_low = np.where(p_range > 0.05)[0][0]
        idx_high = np.where(p_range > 0.8*max(p_range))[0][0]
        
        silent[i] = cdf[idx_low]
        hyper[i] = 1-cdf[idx_high]
        print "silent fraction: ", silent[i]
        print "hyperactive fraction: ", hyper[i]
        
        #silent = np.sum(p_exact[:idx_low])*d_nu
        #hyper = np.sum(p_exact[idx_high:])*d_nu
        print "x value: ", -0.1+0.2/(steps-1)*i
        axsilent.plot(array[i],silent[i],'ko',markersize=3)
        axhyper.plot(array[i],hyper[i],'ro',markersize=3)
        
    
    ax[0].legend()
    ax[0].set_ylabel("1-cdf")
    ax[0].set_xlabel(r"$\displaystyle \nu$ [Hz]")
    
    axsilent.set_xlim([-0.01+array[0],array[-1]+0.01])
    axsilent.set_ylim([0,max(silent)*1.1])
    #axsilent.set_xticks([])
    
    axhyper.set_xlim([-0.01+array[0],array[-1]+0.01])
    axhyper.set_ylim([0,max(hyper)*1.1])
    #axhyper.set_xticks([])

    plt.show(block=False)