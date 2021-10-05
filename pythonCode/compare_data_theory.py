def compare_dat(alpha_0,tau_G):
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
  
  
    d = np.load('../data/firing_rates.npy')
    
    print "mean rate: ", d.mean()
    p_r, p_e = sharkfins(steps=1,rateWnt=[[d.mean()]],alpha_0=[[alpha_0]],tau_G=[[tau_G]],n=[[0]],tau_M=10,mode_calc='exact',mode_stats=3,save=0,file_format='png')

    p_cdf = np.cumsum(p_e)
    p_cdf /= p_cdf[-1]
    
    d.sort()
    dSz = d.size
    print dSz
    d_range = np.zeros(dSz)
    for i in range(dSz):
        d_range[i] = d[i]
        
    d_cdf = np.arange(dSz)/float(dSz)
    
    plt.figure(figsize=(6,4))
    plt.plot(p_r,p_cdf)
    plt.plot(d_range,d_cdf)
    
    #print p_cdf
    #print d_cdf
    
    #plt.figure(figsize=(6,4))
    #plt.hist(d,bins=100,range=[0,10])
    #plt.plot(p_r,p_e*d.size/10.,linewidth=2,color='r')
    plt.xlim([0,10])
    plt.xlabel(r'$\displaystyle \nu$ [Hz]')
    plt.ylabel(r'CDF')
    plt.show(block=False)
    
    
    plt.figure(figsize=(6,4))
    plt.hist(d,bins=100,range=[0,10])
    plt.plot(p_r,p_e*d.size/10.,linewidth=2,color='r')
    plt.xlim([0,10])
    plt.xlabel(r'$\displaystyle \nu$ [Hz]')
    plt.ylabel(r'occurence / $\displaystyle \rho(\nu)$')
    plt.show(block=False)