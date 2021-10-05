


def test_analysis(k=10,j=5,N=None,T=None):
    
    data = {}
    
    if type(N) == int:
        N_sim = N
    else:
        N_sim = 50
        
    if type(T) == int:
        T_sim = T
    else:
        T_sim = 1200
    
    ### scaling with j
    if type(j) == list:
        print "Test scaling with j..."
        
        data['KS'] = np.zeros((k,j[-1]))
        data['KL'] = np.zeros((k,j[-1]))
        
        mean_j = np.zeros((len(j),k))
        var_j = np.zeros((len(j),k))
        
        data_tmp = get_samples_from_theory(alpha_0=0.0,N=N_sim,T=T_sim,n_bin=50,rate=1,prior="mean_field",dftheory=k,dftime=j[-1])
            
            #print data_tmp['KS']
        data['KS'] = data_tmp['KS']
        data['KL'] = data_tmp['KL']
        
        for j_it in range(len(j)):
            #plt.figure()
            #plt.plot(np.arange(k)-0.1,data['KS'][:,:j[j_it]],'ko')
            #plt.plot(np.arange(k)+0.1,data['KL'][:,:j[j_it]],'ro')
            #plt.xlim([-0.5,k-0.5])
            #plt.ylim([0,2])
            #plt.show(block=False)
            
            ### randomize results
            bs_steps=1
            if (bs_steps > 1 and j_it < len(j)):
                mean_tmp = np.zeros((bs_steps,k))
                var_tmp = np.zeros((bs_steps,k))
                for i in range(bs_steps):
                    mean_tmp[i] = np.mean(data['KL'][:,np.random.randint(0,j[-1],j[j_it])],axis=1)
                    var_tmp[i] = np.var(data['KL'][:,np.random.randint(0,j[-1],j[j_it])],axis=1)
                mean_j[j_it] = np.mean(mean_tmp,axis=0)
                var_j[j_it] = np.mean(var_tmp,axis=0)
            else:
                mean_j[j_it] = np.mean(data['KL'][:,:j[j_it]],axis=1)
                var_j[j_it] = np.var(data['KL'][:,:j[j_it]],axis=1)
        
        print "mean: ", mean_j
        print "var: ", var_j
        
        plt.figure()
        ax1 = plt.axes([0.2,0.55,0.7,0.35])
        ax2 = plt.axes([0.2,0.15,0.7,0.35])
        
        ax1.plot(j,mean_j,'k')
        ax2.plot(j,var_j,'r')
        plt.show(block=False)
        
    
    
    ### scaling with N
    if type(N) == list:
        print "Test scaling with N..."
        data['KS'] = np.zeros((len(N),k,j))
        data['KL'] = np.zeros((len(N),k,j))
        
        mean_N = np.zeros((len(N),k))
        var_N = np.zeros((len(N),k))
        
        for N_it in range(len(N)):
            data_tmp = get_samples_from_theory(alpha_0=0.0,N=N[N_it],T=T_sim,n_bin=50,rate=1,prior="mean_field",dftheory=k,dftime=j)
            
            data['KS'][N_it] = data_tmp['KS']
            data['KL'][N_it] = data_tmp['KL']
            
            mean_N[N_it] = np.mean(data['KL'][N_it],axis=1)
            var_N[N_it] = np.var(data['KL'][N_it],axis=1)
            
        plt.figure()
        ax1 = plt.axes([0.2,0.55,0.7,0.35])
        ax2 = plt.axes([0.2,0.15,0.7,0.35])
        
        ax1.plot(N,mean_N,'ko')
        ax2.plot(N,var_N,'ro')
        plt.show(block=False)
        
        plt.figure()
        plt.hist(mean_N[-1])
        plt.hist(data['KL'][-1],bins=np.linspace(0,1,20))
        plt.show(block=False)
                
            
            
            #plt.figure()
            #plt.plot(np.arange(k)-0.1,data['KS'][N_it],'ko')
            #plt.plot(np.arange(k)+0.1,data['KL'][N_it],'ro')
            #plt.xlim([-0.5,k-0.5])
            #plt.ylim([0,2])
            #plt.show(block=False)
        
    return data