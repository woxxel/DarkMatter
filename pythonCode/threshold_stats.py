import imp

imp.load_source('sharks', 'sharkfins.py')
from sharks import *
imp.load_source('single', 'get_samples.py')
from single import *


def threshold_stats(steps=1000,rateWnt=[[1]],alpha_0=[[0,0.02,0.04]],tau_G=[[30]],eps=[[0]],eta=[[0.9]],n=[[0]],Npop=2,drive=0,save=0,file_format='png'):
    
    #steps = 1000
    
    
    #results_stats = sharkfins(steps=steps,rateWnt=[[0,20]],alpha_0=alpha_0,tau_G=[[tau_G]],n=[[0]],drive=drive,mode_calc='exact',mode_stats=1)
    if len(rateWnt[0]) > 1:
        x_ax_string = 'rateWnt'
    if len(eps[0]) > 1:
        x_ax_string = 'eps'
    if len(eta[0]) > 1:
        x_ax_string = 'eta'
    if len(n[0]) > 1:
        x_ax_string = 'n'
    if len(tau_G[0]) > 1:
        x_ax_string = 'tau_G'
    
    two_row = False
    currents = True
    
    x_ax_string = 'n'
    
    results_stats = sharkfins(steps=steps,rateWnt=rateWnt,alpha_0=alpha_0,tau_G=tau_G,n=n,eps=eps,eta=eta,tau_M=10,Npop=Npop,drive=0,mode_calc='exact',mode_stats=1)
    
    print "plotting stuff..."
    
    title_descr = 0
    title_x = -0.2
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    #mpl.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    
    #print results_stats['nu_c_idx']
    #print results_stats['nu_c']
    
    
    #if not x_ax_string == 'n':
    if two_row:
        fig, ax = plt.subplots(2,3,figsize=(7.5,3.5))
    else:
        fig, ax = plt.subplots(2,3,figsize=(7.5,2))
        fig.delaxes(ax[1,0])
        fig.delaxes(ax[1,1])
        fig.delaxes(ax[1,2])
        ax[0,0].set_position([0.1,0.25,0.22,0.6])
        ax[0,1].set_position([0.415,0.25,0.22,0.6])
        ax[0,2].set_position([0.76,0.25,0.22,0.6])
        
        
    
    #big_ax = fig.add_subplot(111)
    big_ax = plt.axes([0.1,0.1,0.8,0.8])
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])
    
    if x_ax_string == 'rateWnt':
        big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    elif x_ax_string == 'eps':
        big_ax.set_xlabel(r'$\displaystyle \varepsilon$')
    elif x_ax_string == 'n':
        big_ax.set_xlabel(r'$\displaystyle b$')
    
    #print results_stats['nu_implausible']
    #print results_stats
    nu_imp_idx = np.zeros(len(results_stats['gamma']))
    
    #print results_stats['gamma
    for a in range(len(results_stats['gamma'])):
        #print results_stats[x_ax_string]
        
        try:
            nu_imp_idx[a] = np.where(results_stats[x_ax_string][a]==results_stats['nu_implausible'][a])[0][0]
        except:
            nu_imp_idx[a] = -1
        #print "implausible: ", nu_imp_idx[a], results_stats[x_ax_string][a][nu_imp_idx[a]]
        #if not x_ax_string == 'n':
        if two_row:
            col = a/float(len(results_stats['gamma']))
            
            #if x_ax_string == 'rateWnt':
                #ax[1,1].plot(results_stats[x_ax_string][a],results_stats['gamma'][a]**2,color=(col,0,0),linestyle=':')
                #ax[1,1].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]],results_stats['gamma'][a][:nu_imp_idx[a]]**2,color=(col,0,0))
            #if x_ax_string in ['eps','eta','n','tau_G']:
            ax[1,1].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]],results_stats['gamma'][a][:nu_imp_idx[a]]**2,color='r')
            if Npop == 2:
                ax[1,1].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]],results_stats['gamma_exc'][a][:nu_imp_idx[a]]**2,color='k')
            
            ax[1,0].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi'][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color='r')
            if Npop == 2:
                ax[1,0].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi_exc'][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color='k')
                
            #mask = np.where(results_stats['chi'][a] > 0)
            
            #print "nu_c_idx: ", results_stats['nu_c_idx'][a]
            #print results_stats['chi'][a]
            #print results_stats['chi'][a]
            #print results_stats['chi'][a][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]]
            #if not (results_stats['nu_c_idx'][a][0] == results_stats['nu_c_idx'][a][1]):
                
                #if two_row:
                    #ax[1,0].plot(results_stats[x_ax_string][a][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi'][a][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color=(col,0,0),linestyle=':')
                    #ax[1,0].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi'][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color=(col,0,0))
                #else:
                ##elif x_ax_string in ['eps','eta','n','tau_G']:
                    #ax[1,0].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi'][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color='r')
                    #ax[1,0].plot(results_stats[x_ax_string][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],results_stats['chi_exc'][a][:nu_imp_idx[a]][results_stats['nu_c_idx'][a][0]:results_stats['nu_c_idx'][a][1]],color='k')
    
    
    
    #if x_ax_string == 'rateWnt':
    if not currents:
        ax[0,0].plot(results_stats[x_ax_string][0],results_stats['q'][0],'k:')
        ax[0,0].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['q'][0][:nu_imp_idx[0]],'k',label=r'$\displaystyle \alpha_0 = 0$')
        if (len(results_stats[x_ax_string]) > 1):
            ax[0,0].plot(results_stats[x_ax_string][2],results_stats['q'][2],'r:')
            ax[0,0].plot(results_stats[x_ax_string][2][:nu_imp_idx[2]],results_stats['q'][2][:nu_imp_idx[2]],'r',label=r'$\displaystyle \alpha_0 = 0.04$')
        
        ax[0,0].set_ylabel(r'$\displaystyle q\,$[Hz$\displaystyle ^2$]')
    
        ax[0,0].plot(results_stats[x_ax_string][0,:0.1*steps],results_stats[x_ax_string][0,:0.1*steps]*results_stats['rate_max']/np.sqrt(2),'r--',linewidth=0.5)
        if title_descr:
            ax[0,0].set_title('a) second moment',position=(0.1,1.05))#,loc='left'
        else:
            ax[0,0].set_title(r'a)',position=(title_x,1.05))
        
        
        ax[0,0].plot(results_stats[x_ax_string][0],results_stats[x_ax_string][0]**2,'k--',linewidth=0.5)
        ax[0,0].set_xlim([0,15])
        ax[0,0].set_ylim([0,200])
        plt.setp(ax[0,0],yticks=np.linspace(0,200,5))
        ax[0,0].text(8,43,r'$\displaystyle q\approx\bar{\nu}^2$',fontsize=10)
        ax[0,0].legend(prop={'size':10},bbox_to_anchor=(0,1.2),loc='lower right')
        
    else:
    #elif x_ax_string in ['eps','eta','n','tau_G']:
        #elif x_ax_string in ['eps','eta','n']:
        #ax[0,0].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['q'][0][:nu_imp_idx[0]],'r',label='inh.')
        #ax[0,0].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['q_exc'][0][:nu_imp_idx[0]],'k',label='exc.')#r'$\displaystyle \alpha_0 = 0$')
        
        #ax[0].legend(prop={'size':10},loc='lower left')
        
        ax[0,0].legend(prop={'size':10},loc='lower left')
        #ax[0,0].set_ylim([2,4])
        #plt.setp(ax[0,0],yticks=np.linspace(2,4,5))

        tau_M = 0.010
        
        tau_A = 0.005
        #tau_G = 0.030
        tau_N = 0.200
        
        J_EE = results_stats['eta'] * results_stats['eps']*tau_M
        J_EI = np.sqrt(1-(results_stats['eta']*results_stats['eps'])**2)*tau_M

        
        var_V_A = J_EE**2 * results_stats['rateWnt'][0] / (tau_A + tau_M) * ( (1-results_stats['n'][0])**2/2 + (1-results_stats['n'][0])*results_stats['n'][0]*tau_A / (tau_A + tau_N) )
        #from excitatory NMDA synapses
        var_V_N = J_EE**2 * results_stats['rateWnt'][0] / (tau_N + tau_M) * (results_stats['n'][0]**2/2 + (1-results_stats['n'][0])*results_stats['n'][0]*tau_N / (tau_A + tau_N) );
        #from inhibitory GABA synapses
        var_V_G = J_EI**2 * results_stats['rateWnt'][0] * 0.5 / (results_stats['tau_G'][0]/1000. + tau_M)
        
        #return results_stats, var_V_A, var_V_G
        #print results_stats['n']
        #print var_V_N
        #print var_V_G
        ax[0,0].plot(results_stats[x_ax_string][0],np.sqrt(var_V_G[0]),'k--',label=r'$\displaystyle \sigma_{V_E^G}$ (GABA)')
        ax[0,0].plot(results_stats[x_ax_string][0],np.sqrt(var_V_A[0]),'k:',label=r'$\displaystyle \sigma_{V_E^A}$ (AMPA)')
        #ax[0,0].plot(results_stats[x_ax_string][0],np.sqrt(var_V_G[0]+var_V_A[0]),'k-')
        if any(max(n)) > 0:
            ax[0,0].plot(results_stats[x_ax_string][0],np.sqrt(var_V_N[0]),'k-.',label=r'$\displaystyle \sigma_{V_E^N}$ (NMDA)')
            ax[0,0].legend(prop={'size':10},bbox_to_anchor=(0.05,1.2),loc='upper left',ncol=1)
        else:
            ax[0,0].text(0.2,0.005,r'$\displaystyle \sigma_{V_{EE}}$ (AMPA) $\displaystyle \hat{=} \sigma_{V_E^A}$',fontsize=10)
            ax[0,0].text(0.1,0.038,r'$\displaystyle \sigma_{V_{EI}}$ (GABA) $\displaystyle \hat{=} \sigma_{V_E^G}$',fontsize=10)
        
        if x_ax_string == 'rateWnt':
            plt.setp(ax[0,0],xticks=np.linspace(0,10,6),yticks=np.linspace(0,0.15,4),xlim=[0,results_stats[x_ax_string][0][-1]],ylim=[0,0.15])
        else:
            #plt.setp(ax[0,0],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,4),xlim=[0,results_stats[x_ax_string][0][-1]])#,ylim=[0,0.075])
            plt.setp(ax[0,0],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,6),xlim=[0,results_stats[x_ax_string][0][-1]])
        ax[0,0].set_ylabel(r'$\displaystyle I$')
        ax[0,0].set_title(r'a)',position=(title_x,1.05))
        
    
    
    
    
    #if x_ax_string == 'rateWnt':
    if not currents:
        ax[0,1].plot(results_stats[x_ax_string][0],results_stats['q'][0],'k')
        if (len(results_stats[x_ax_string]) > 1):
            ax[0,1].plot(results_stats[x_ax_string][2],results_stats['q'][2],'r')
        
        ax[0,1].set_ylabel(r'$\displaystyle q\,$[Hz$\displaystyle ^2$]')
        if title_descr:
            ax[0,1].set_title(r'b) the low $\displaystyle \bar{\nu}$ limit',position=(0.1,1.05))#loc='left')
        else:
            ax[0,1].set_title(r'b)',position=(title_x,1.05))
        
        ax[0,1].set_xlim([0,1])
        ax[0,1].set_ylim([0,5])
        ax[0,1].text(0.1,3.5,r'$\displaystyle q\approx \frac{\bar{\nu}\nu_{max}}{\sqrt{2}}$',fontsize=10)
        
        ax[0,1].plot(results_stats[x_ax_string][0,:0.012*steps],results_stats[x_ax_string][0,:0.012*steps]*results_stats['rate_max']/np.sqrt(2),'r--',linewidth=0.5)
        ax[0,1].plot(results_stats[x_ax_string][0],results_stats[x_ax_string][0]**2,'k--',linewidth=0.5)
        #plt.setp(ax[0,1],xticks=np.linspace(0,1,6))
        plt.setp(ax[0,1],xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
    else:
    #elif x_ax_string in ['eps','eta','n','tau_G']:
        ax[0,1].plot(results_stats[x_ax_string][0],results_stats['alpha'][0],'r:')
        if (len(results_stats[x_ax_string]) > 1):
            ax[0,1].plot(results_stats[x_ax_string][2],results_stats['alpha'][2],'r:')#,label=r'$\displaystyle \alpha$')
        ax[0,1].plot(results_stats[x_ax_string][0][nu_imp_idx[0]:],results_stats['sigma_V'][0][nu_imp_idx[0]:],'r:')
        
        ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha'][0][:nu_imp_idx[0]],'r-',label='inh.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
        if Npop == 2:
            ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha_exc'][0][:nu_imp_idx[0]],'k-',label='exc.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
        
        #if (len(results_stats[x_ax_string]) > 1):
            #ax[0,1].plot(results_stats[x_ax_string][2][:nu_imp_idx[2]],results_stats['alpha'][2][:nu_imp_idx[2]],'r-')#,label=r'$\displaystyle \alpha$')
        ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V'][0][:nu_imp_idx[0]],'r--')#,label=r'$\displaystyle \sigma_V$')
        if Npop == 2:
            ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V_exc'][0][:nu_imp_idx[0]],'k--')#,label=r'$\displaystyle \sigma_V$')
        
        #ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['I_balance'][0][:nu_imp_idx[0]],'-',color=[0.7,0.7,0.7],label=r'$\displaystyle I_{balance}$')

        if (len(results_stats[x_ax_string]) > 1):
            if not np.isnan(results_stats['nu_c'][2]):
                ax[0,2].plot([results_stats['nu_c'][2],results_stats['nu_c'][2]],[0,1],'k:')
        
        
        ax[0,1].text(0.1,0.015,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)
        if x_ax_string == 'n':
            ax[0,1].text(0.05,0.05,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k^A}^2 + \sigma_{V_k^N}^2 + \sigma_{V_k^G}^2}$',fontsize=10)
        else:
            ax[0,1].text(0.05,0.043,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k^A}^2 + \sigma_{V_k^G}^2}$',fontsize=10)

        #ax[0,1].text(0.05,0.016,r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$',fontsize=10)
        if (len(results_stats[x_ax_string]) > 1):
            if not np.isnan(results_stats['nu_c'][2]):
                ax[0,1].text(results_stats['nu_c'][2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)
        
        
        #ax[0,1].set_ylim([0.0,0.05])
        #plt.setp(ax[0,1],yticks=np.linspace(0.02,0.05,4),xticks=np.linspace(0,0.8,5))
        ax[0,1].set_ylabel(r'$\displaystyle I$')
        
        if x_ax_string == 'rateWnt':
            plt.setp(ax[0,1],xticks=np.linspace(0,10,6),yticks=np.linspace(0,0.15,4),xlim=[0,results_stats[x_ax_string][0][-1]],ylim=[0,0.15])
        else:
            #plt.setp(ax[0,1],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,4),xlim=[0,results_stats[x_ax_string][0][-1]])
            plt.setp(ax[0,1],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,6),xlim=[0,results_stats[x_ax_string][0][-1]])
        #print "lim: ", results_stats[x_ax_string][0][-1]
        #ax[0,1].set_xlim([0,results_stats[x_ax_string][0][-1]])
        #ax[0,1].legend(prop={'size':10},loc='lower right')
        #ax[0,2].set_yticks(np.linspace(0,0.15,4))
        #ax[0,2].set_yticklabels(np.linspace(0,0.15,4),fontsize=12)
        if title_descr:
            ax[0,1].set_title('b) variances',position=(0.1,1.05))#,loc='left')
        else:
            ax[0,1].set_title(r'b)',position=(title_x,1.05))
    
    
    
    
    
    #ax[0,2].plot(results_stats[x_ax_string][0],np.zeros(len(results_stats[x_ax_string][0])),'k:',linewidth=0.5)
    
    #if x_ax_string == 'rateWnt':
        
            
    #if x_ax_string in ['eps','eta','n','tau_G']:
    if currents:
        ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],-results_stats['I_balance'][0][:nu_imp_idx[0]],'-',color='r',label='inh.')
        if Npop == 2:
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],-results_stats['I_balance_exc'][0][:nu_imp_idx[0]],'-',color='k',label='exc.')
        
        if x_ax_string == 'rateWnt':
            plt.setp(ax[0,2],xticks=np.linspace(0,10,6),yticks=np.linspace(-0.15,0,4),xlim=[0,results_stats[x_ax_string][0][-1]],ylim=[-0.15,0])
        else:
            plt.setp(ax[0,2],xticks=np.linspace(0,1,6),yticks=np.linspace(-0.15,0,4),xlim=[0,results_stats[x_ax_string][0][-1]])
        ax[0,2].legend(prop={'size':10},bbox_to_anchor=(0.95,1.1),loc='upper right')
        

    
    
    ax[0,2].text(6,-0.23,r'$\displaystyle \bar{I}_0-\Psi_0$',fontsize=10)
    
    
    try_square = False
    #if x_ax_string == 'rateWnt':
    if not currents:
        
        if try_square:
            ax[0,2].plot(results_stats[x_ax_string][0],0.01**2*results_stats[x_ax_string][0]**2,'r-')
            
            ax[0,2].plot(results_stats[x_ax_string][0][nu_imp_idx[0]:],results_stats['sigma_V'][0][nu_imp_idx[0]:]**2,'k:')
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V'][0][:nu_imp_idx[0]]**2,'k--')
            
            ax[0,2].plot(results_stats[x_ax_string][0][nu_imp_idx[0]:],results_stats['alpha'][0][nu_imp_idx[0]:]**2,'k:',label=r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha'][0][:nu_imp_idx[0]]**2,'k-')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
            ax[0,2].set_ylim([0,max(results_stats['sigma_V'][0])**2])
        
        else:
                
            ax[0,2].plot(results_stats[x_ax_string][0],-results_stats['I_balance'][0],':',color=[0.7,0.7,0.7])
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],-results_stats['I_balance'][0][:nu_imp_idx[0]],'-',color=[0.7,0.7,0.7],label=r'$\displaystyle I_{balance}$')
            
            if (len(results_stats[x_ax_string]) > 1):
                if not np.isnan(results_stats['nu_c'][2]):
                    ax[0,2].plot([results_stats['nu_c'][2],results_stats['nu_c'][2]],[0,1],'k:')
            if (len(results_stats[x_ax_string]) > 1):
                if not np.isnan(results_stats['nu_c'][2]):
                    ax[0,2].text(results_stats['nu_c'][2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)
                    
            ax[0,2].plot(results_stats[x_ax_string][0][nu_imp_idx[0]:],results_stats['sigma_V'][0][nu_imp_idx[0]:],'k:')
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V'][0][:nu_imp_idx[0]],'k--',label=r'$\displaystyle \sigma_V$')
            
            ax[0,2].plot(results_stats[x_ax_string][0],results_stats['alpha'][0],'k:')
            ax[0,2].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha'][0][:nu_imp_idx[0]],'k-',label=r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
            
            if (len(results_stats[x_ax_string]) > 1):
                ax[0,2].plot(results_stats[x_ax_string][2],results_stats['alpha'][2],'r:')#,label=r'$\displaystyle \alpha$')
                ax[0,2].plot(results_stats[x_ax_string][2][:nu_imp_idx[2]],results_stats['alpha'][2][:nu_imp_idx[2]],'r-')#,label=r'$\displaystyle \alpha$')
            
            
            ax[0,2].text(6,0.17,r'$\displaystyle \sigma_{V_k}$',fontsize=10)
            ax[0,2].text(8,0.02,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)
        
        #ax[0,2].set_xlim([0,15])
        plt.setp(ax[0,2],xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
        #ax[0,2].set_ylim([-0.3,0.25])
    #elif x_ax_string in ['eps','eta','n']:
        #ax[0,2].set_xlim([0,0.8])
        #ax[0,2].set_ylim([-0.15,-0.05])
        #plt.setp(ax[0,2],yticks=np.linspace(-0.15,-0.05,3),xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
                 
    ax[0,2].set_ylabel(r'$\displaystyle \bar{I} - \Psi_0$')
    #plt.setp(ax[0,2],yticks=np.linspace(0,0.15,4))
    #ax[0,2].legend(prop={'size':10},loc=2)
    #ax[0,2].set_yticks(np.linspace(0,0.15,4))
    #ax[0,2].set_yticklabels(np.linspace(0,0.15,4),fontsize=12)
    if title_descr:
        ax[0,2].set_title('c) variances',position=(0.1,1.05))#,loc='left')
    else:
        ax[0,2].set_title(r'c)',position=(title_x,1.05))
    
    
    #if not x_ax_string == 'n':
    if two_row:
        
        ax[1,1].plot(results_stats[x_ax_string][0],np.ones(len(results_stats[x_ax_string][0])),'k:',linewidth=0.5)
        #ax[1,1].plot(results_stats[x_ax_string][0],results_stats['gamma'][0]**2-1,'k')
        #ax[1,1].plot(results_stats[x_ax_string][1],results_stats['gamma'][1]**2-1,color=[0.5,0,0])
        #ax[1,1].plot(results_stats[x_ax_string][2],results_stats['gamma'][2]**2-1,color=[1,0,0]])
        if title_descr:
            ax[1,1].set_title(r'e)$\displaystyle\quad$ DM exponent $\displaystyle \gamma$',position=(0.1,1.05))#,loc='left')
        else:
            ax[1,1].set_title(r'e)',position=(title_x,1.05))
                
        #ax[1,1].plot()
        if x_ax_string == 'rateWnt':
            ax[1,1].set_xlim([0,15])
            ax[1,1].set_ylim([0,7])
            
            
                
        elif x_ax_string in ['eps','eta']:
            plt.setp(ax[1,1],xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
            ax[1,1].set_title(r'e)',position=(title_x,1.05))
            
        ax[1,1].set_ylabel(r'$\displaystyle \gamma^2$')
        
        
        
        
        
        
        #ax[1,0].plot(results_stats[x_ax_string][0][results_stats['nu_c_idx'][0]:],results_stats['chi'][0][results_stats['nu_c_idx'][a]:],'k')
        #ax[1,0].plot(results_stats[x_ax_string][1][results_stats['nu_c_idx'][1]:],results_stats['chi'][a][results_stats['nu_c_idx'][a]:],'k')
        #ax[1,0].plot(results_stats[x_ax_string][2][results_stats['nu_c_idx'][2]:],results_stats['chi'][a][results_stats['nu_c_idx'][a]:],'k')
        
        #if x_ax_string == 'rateWnt':
        #if two_row:
            
        #ax[1,0].text(results_stats[x_ax_string][0,0.1*steps],0.5*results_stats['chi'][0,0.1*steps],r'$\displaystyle \alpha_0 = 0$',fontsize=8)
        #if (len(results_stats[x_ax_string]) > 1):
            #ax[1,0].text(results_stats[x_ax_string][1,0.07*steps],1.1*results_stats['chi'][1,0.07*steps],r'$\displaystyle \alpha_0 = 0.02$',fontsize=8)
            #ax[1,0].text(results_stats[x_ax_string][2,0.3*steps],results_stats['chi'][2,0.3*steps],r'$\displaystyle \alpha_0 = 0.04$',fontsize=8)
        
        #ax[1,0].set_xlim([0,15])
        #if (max(results_stats['chi']) > 10)
        ax[1,0].set_ylim([0,10])
        #elif x_ax_string in ['eps','eta','n','tau_G']:
            #plt.setp(ax[1,0],xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
            
        ax[1,0].set_ylabel(r'$\displaystyle \chi$')
        plt.setp(ax[1,0],xticks=np.linspace(0,0.8,5),yticks=np.linspace(0,10,3),xlim=[0,results_stats[x_ax_string][0][-1]])
        if title_descr:
            ax[1,0].set_title(r'd) skewness coeff. $\displaystyle \chi$',position=(0.125,1.05))#,loc='left')#,fontsize=12)
        else:
            ax[1,0].set_title(r'd)',position=(title_x,1.05))
        
        
        
        
        #if x_ax_string == 'rateWnt' and not try_square:
            #results_bounds = sharkfins(steps=100,rateWnt=[[0,20]],alpha_0=[[0,0.2]],tau_G=[[tau_G]],n=[[0]],drive=drive,mode_calc='exact',mode_stats=2)
            
            #mask_DM = ~np.isnan(results_bounds['DM_trans'])
            #mask_no_peak = ~np.isnan(results_bounds['no_peak_trans'])
            #mask_inc = ~np.isnan(results_bounds['inc_trans'])
            
            #ax[1,2].plot(results_bounds['DM_trans'][mask_DM],results_bounds['alpha_0'][mask_DM],'r-',label=r'$\displaystyle \bar{\nu}_{DM}$')
            #ax[1,2].plot(results_bounds['no_peak_trans'][mask_no_peak],results_bounds['alpha_0'][mask_no_peak],'k--',label=r'$\displaystyle \bar{\nu}_{no\,peak}$')
            #ax[1,2].plot(results_bounds['inc_trans'][mask_inc],results_bounds['alpha_0'][mask_inc],'k',label=r'$\displaystyle \bar{\nu}_{inc}$')
            #ax[1,2].plot(results_bounds['nu_implausible'][mask_inc],results_bounds['alpha_0'][mask_inc],'k:',label=r'$\displaystyle \bar{\nu}_{imp}$')
            
            #if x_ax_string == 'rateWnt':
                #ax[1,2].set_xlim([0,20])
                #ax[1,2].set_ylim([0,0.2])
            #elif x_ax_string in ['eps','eta','n']:
                #plt.setp(ax[1,2],xticks=np.linspace(0,0.8,5),xlim=[0,results_stats[x_ax_string][0][-1]])
            
            #ax[1,2].set_yticks(np.linspace(0,0.2,6))
            #ax[1,2].set_ylabel(r'$\displaystyle \alpha_0$')
            
            #if title_descr:
                #ax[1,2].set_title(r'f) boundaries',position=(0.1,1.05))#,loc='left')
            #else:
                #ax[1,2].set_title(r'f)',position=(title_x,1.05))
            #ax[1,2].legend(prop={'size':10},ncol=2,bbox_to_anchor=(-0.05,1.3),loc='upper left')
            
        #if x_ax_string in ['eps','eta','n','tau_G']:
            ##print "eps: ", results_stats['eps'][0]
            ##print "eta: ", results_stats['eta'][0][0]
            #I_I_per_nu = np.sqrt(1-results_stats['eps'][0]**2) - results_stats['eps'][0]
            
            #eta = [0.9,0.6,0.2]
            #for i in range(len(eta)):
                #I_E_per_nu = np.sqrt(1-(eta[i]*results_stats['eps'][0])**2) - eta[i]*results_stats['eps'][0]
                #col = float(i)/len(eta)
                #ax[1,2].plot(results_stats[x_ax_string][0],I_E_per_nu,color=(col,col,col),label=r'$\displaystyle \eta = %3.1g$'%eta[i])
            
            #ax[1,2].plot(results_stats[x_ax_string][0],I_I_per_nu,'r')
            #ax[1,2].legend(prop={'size':10},loc='lower left')
            #ax[1,2].set_ylabel(r'$\displaystyle I^{ext} / \bar{\nu}$')
            
            #plt.setp(ax[1,2],xticks=np.linspace(0,0.8,5),yticks=np.linspace(0,1,3),xlim=[0,results_stats[x_ax_string][0][-1]])
            
            #ax[1,2].set_title(r'f)',position=(title_x,1.05))
        
        plt.subplots_adjust(left=0.075, bottom=0.11, right=0.99, top=0.925, wspace=0.4, hspace=0.35)
    mpl.rcParams['axes.titlesize'] = 12
    #plt.show()
    
    ## for the last plot: need new simulation, that samples nu-alpha_0 phase space for boundaries with lower sampling rate
    
    
    for i in range(2):
        for j in range(3):
            ax[i,j].spines['right'].set_color('none')
            ax[i,j].yaxis.set_ticks_position('left')
            ax[i,j].spines['top'].set_color('none')
            ax[i,j].xaxis.set_ticks_position('bottom')
        
        #ax[i,j].tick_params(axis='both', which='major', labelsize=16)
        #ax[i,j].tick_params(axis='both', which='minor', labelsize=16)
        #ax[i,j].set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    #mpl.rcParams['xtick.labelsize'] = 12
    #plt.rcParams['xtick.labelsize'] = 'x-large'
    
    if save:
        if Npop == 1:
            sv_name = './../paper draft/inhib only/pics/threshold_stats_%s_drive_%d.%s' % (x_ax_string,drive,file_format)
        if Npop == 2:
            sv_name = './../paper draft/two_pop/pics/threshold_stats_%s_drive_%d.%s' % (x_ax_string,drive,file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print 'Figure saved as "%s"' % sv_name
    
    plt.show(block=False)
    
    #print results_stats.keys()
    #y = results_stats['sigma_V'][0]/results_stats['I_balance'][0]
    #print results_stats['sigma_V'][0]
    #print results_stats['I_balance'][0]
    #print y
    #plt.figure()
    #plt.plot(results_stats[x_ax_string][0],y)
    #plt.plot(results_stats[x_ax_string][0],results_stats['alpha_raw'][0]/max(results_stats['alpha'][0]),'r--')
    #plt.plot(results_stats[x_ax_string][0],results_stats['alpha'][0]/max(results_stats['alpha'][0]),'r-')
    
    #plt.show(block=False)


def plot_distr(alpha_0=[0,0.04],rateWnt=[0.5,2,10],tau_M=0.010,save=0,file_format='png'):
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    
    title_x = -0.1
    title_y = 1.05
    
    plt.figure(figsize=(7.5,4))
    
    nu_border = [5,15,25]
    
    v_bottom = 0.1
    v_top = 0.925
    
    v_spaces = 0.1
    v_spaces_small = 0.02
    v_plots = ((v_top-v_bottom)-2*v_spaces)/3.
    print "v plots: ", v_plots
    
    big_ax1 = plt.axes([0.1,0.05,0.35,0.9])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax1.spines['top'].set_visible(False)
    big_ax1.spines['right'].set_visible(False)
    big_ax1.spines['bottom'].set_visible(False)
    big_ax1.spines['left'].set_visible(False)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_xlabel(r'$\displaystyle I$')
    
    big_ax1 = plt.axes([0.05,0.05,0.35,0.75 ])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax1.spines['top'].set_visible(False)
    big_ax1.spines['right'].set_visible(False)
    big_ax1.spines['bottom'].set_visible(False)
    big_ax1.spines['left'].set_visible(False)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_ylabel(r'$\displaystyle \rho(I)$',color='r')
    
    big_ax1 = plt.axes([0.05,0.25,0.35,0.75])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax1.spines['top'].set_visible(False)
    big_ax1.spines['right'].set_visible(False)
    big_ax1.spines['bottom'].set_visible(False)
    big_ax1.spines['left'].set_visible(False)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_ylabel(r'$\displaystyle \nu$ [Hz]')
    
    big_ax2 = plt.axes([0.6,0.075,0.35,0.9])
    big_ax2.set_axis_bgcolor('none')
    big_ax2.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax2.spines['top'].set_visible(False)
    big_ax2.spines['right'].set_visible(False)
    big_ax2.spines['bottom'].set_visible(False)
    big_ax2.spines['left'].set_visible(False)
    plt.setp(big_ax2,xticks=[],yticks=[])
    big_ax2.set_xlabel(r'$\displaystyle \rho(\nu)$')
    
    for i in range(len(rateWnt)):
        for j in range(len(alpha_0)):
            results = get_samples_from_theory(tau_M=tau_M,T=1000,dftheory=0,dftime=0,p_theory=2,plot=1,rate=rateWnt[i],alpha_0=alpha_0[j])
            #print results
            #ax1 = plt.subplot(gs[2*i,0])
            if (j == 0):
                border_y = (len(rateWnt)-i-1)*v_plots + (len(rateWnt)-i-1)*v_spaces + v_bottom
                #ax1 = plt.axes([0.1,border_y+(v_plots+v_spaces_small)/2.,0.33,(v_plots-v_spaces_small)/2.])
                #ax2 = plt.axes([0.1,border_y,0.33,(v_plots-v_spaces_small)/2.])
                ax2 = plt.axes([0.1,border_y,0.45,v_plots])
                ax1 = ax2.twinx()
                ax3 = plt.axes([0.6,border_y,0.15,v_plots])
                
                ax2.set_xticks(np.linspace(-0.4,0.0,5))
                ax2.set_yticks([])
                ##ax2.set_xlabel(r'$\displaystyle I$')
                #ax2.set_ylabel(r'$\displaystyle \rho(I)$')
                ax2.spines['right'].set_color('none')
                ax2.yaxis.set_ticks_position('left')
                ax2.spines['top'].set_color('none')
                ax2.xaxis.set_ticks_position('bottom')
                
                #ax1.set_xticks([])
                #ax1.set_ylabel(r'$\displaystyle \nu(I)\,$[Hz]')
                ax1.spines['right'].set_color('none')
                ax1.yaxis.set_ticks_position('left')
                ax1.spines['top'].set_color('none')
                ax1.xaxis.set_ticks_position('bottom')

                ax3.set_xticks([])
                ax3.set_yticks([])
                #ax3.set_xlabel(r'$\displaystyle \nu\,$[Hz]')
                #ax3.set_ylabel(r'$\displaystyle \rho(\nu)$')
                ax3.spines['right'].set_color('none')
                ax3.yaxis.set_ticks_position('left')
                ax3.spines['top'].set_color('none')
                ax3.xaxis.set_ticks_position('bottom')
                
                #factor = 1-1./math.pi
                factor = 0.6
                
                nu_max = max(results['f_I'])
                print nu_max
                ax1.plot([-results['sigma_V'],0],[factor*nu_max,factor*nu_max],'k',linewidth=2)#,label=r'$\displaystyle \sigma_V$')
                ax1.plot(results['I_range'],results['f_I'],'k')
                ax1.plot([0,0],[0,nu_max],'k--',linewidth=0.5)
                
                ax1.set_yticks(np.linspace(0,25,6))
                
                
                
                #ax1.legend(prop={'size':10},loc=3)
                
                I_max = max(results['I_distr'])
                I_tmp = np.copy(I_max)
                ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r',linewidth=2,label=r'$\displaystyle \alpha_I$')
                ax2.plot(results['I_range'],results['I_distr'],'r')
                ax2.plot([results['I'],results['I']],[0,I_tmp],'r--',linewidth=0.5)
                #ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r',linewidth=2,label=r'$\displaystyle \alpha_I$')
                #ax2.plot(results['I_range'],results['I_distr'],'k')
                #ax2.plot([results['I'],results['I']],[0,I_max],'k--',linewidth=0.5)
                #ax2.plot([0,0],[0,1.1*I_max],'k--',linewidth=0.5)
                
                
                
                
                #ax2.legend(prop={'size':10},loc=4)
                
                #print results['p_range']
                rho_at_mean = results['p_exact'][np.argmin(abs(results['p_range']-rateWnt[i]))]
                rho_max1 = max(results['p_exact'])
                
                #ax3.plot([rateWnt[i],rateWnt[i]],[0,rho_at_mean],'k--')
                #ax3.plot(results['p_range'],results['p_exact'],'k',label='exact solution')
                #ax3.plot(results['p_range'],results['p_approx'],'r--',linewidth=1,label='approx. solution')
                #ax3.set_xlim([0,nu_border[i]])
                ax3.plot([0,rho_at_mean],[rateWnt[i],rateWnt[i]],'k--')
                ax3.plot(results['p_exact'],results['p_range'],'k',label='exact solution')
                ax3.plot(results['p_approx'],results['p_range'],'r--',linewidth=1,label='approx. solution')
                
                chi1 = np.copy(results['chi'])
                #print chi1
                if (i==0):
                    ax1.set_title(r'a)',position=(title_x,title_y))#,loc='left')
                    ax3.set_title(r'b)',position=(title_x,title_y))# homog.: $\displaystyle \alpha_0 = 0$',loc='left')
                    ax1.text(0.01,1,r'$\displaystyle \Psi$',fontsize=12)
                    y_border = 12
                else:
                    y_border = 25
                    
                if (i==1):
                    #if results['sigma_V'] > 0.05:
                    ax1.text(-results['sigma_V']/2.-0.01,factor*nu_max*0.7,r'$\displaystyle \sigma_V$',fontsize=10)
                    #if results['alpha'] > 0.05:
                    ax2.text(results['I']-results['alpha']/2.-0.01,factor*I_max*0.7,r'$\displaystyle \alpha_I$',fontsize=10)
                
                if (i==2):
                    ax1.text(results['I']+0.01,1,r'$\displaystyle \bar{I}_0$',fontsize=10)
                    
                ax1.set_ylim([0,y_border])
                ax2.set_ylim([0,1.1*I_max*y_border/25.])
                ax3.set_ylim([0,y_border])
                
                ax1.text(-0.45,0.9*y_border,r'$\displaystyle \bar{\nu}=%g\,$Hz'%rateWnt[i],fontsize=12)
                
                
                    #ax3.text(rateWnt[i]+0.2,rho_at_mean,r'$\displaystyle \bar{\nu}$',fontsize=12)
                #else:
                #ax3.text(rho_at_mean,rateWnt[i]+0.05*nu_border[i],r'$\displaystyle \bar{\nu}$',fontsize=12)
            if (j == 1):
                #ax1.plot(results['I_range'],results['f_I'],'k')
                I_max = max(results['I_distr'])
                #ax2.plot(results['I_range'],results['I_distr']*I_tmp/I_max,'k:',linewidth=0.5)
                ax2.plot(results['I_range'],results['I_distr']*I_tmp/I_max,'r:',linewidth=0.5)
                ax2.plot([results['I'],results['I']],[0,I_tmp],'r:',linewidth=0.5)
                #factor = (1-1./np.exp(1))*0.9
                #ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r--',linewidth=2)
                ax4 = plt.axes([0.8,border_y,0.15,v_plots])
                
                ax4.set_xticks([])
                ax4.set_yticks([])
                #ax4.set_xlabel(r'$\displaystyle \nu\,$[Hz]')
                #ax4.set_ylabel(r'$\displaystyle \rho(\nu)$')
                ax4.spines['right'].set_color('none')
                ax4.yaxis.set_ticks_position('left')
                ax4.spines['top'].set_color('none')
                ax4.xaxis.set_ticks_position('bottom')
                
                rho_at_mean = results['p_approx'][np.argmin(abs(results['p_range']-rateWnt[i]))]
                rho_max2 = max(results['p_exact'])
                #print results['p_exact']
                if (rho_max2 == results['p_exact'][1]):
                    rho_max2 /= 3
                rho_max = max(rho_max1,rho_max2)
                #ax4.plot([rateWnt[i],rateWnt[i]],[0,rho_at_mean],'k--')
                #ax4.plot(results['p_range'],results['p_exact'],'k')
                #ax4.plot(results['p_range'],results['p_approx'],'r--')
                #ax4.set_xlim([0,nu_border[i]])
                ax4.plot([0,rho_at_mean],[rateWnt[i],rateWnt[i]],'k--')
                ax4.plot(results['p_exact'],results['p_range'],'k')
                ax4.plot(results['p_approx'],results['p_range'],'r--')
                
                if (i==0):
                    ax4.set_title(r'c)',position=(title_x,title_y))# inhom.: $\displaystyle \alpha_0 = 0.04$',loc='left')
                    ax4.legend(prop={'size':10},bbox_to_anchor=(-0.4,0.85),loc='lower left')
                
                ax4.set_ylim([0,y_border])
                
                ax3.text(0.6*rho_max,y_border*3/4.,r'$\displaystyle \chi \approx %4.2f$'%chi1,bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=10)
                ax4.text(0.6*rho_max,y_border*3/4.,r'$\displaystyle \chi \approx %4.2f$'%results['chi'],bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=10)
                
                #else:
                    #ax3.text(12.5,0.7*rho_max,r'$\displaystyle \chi \approx %4.2g$'%chi1,bbox={'facecolor':'white','alpha':0.9,'pad':5})
                    #ax4.text(12.5,0.7*rho_max,r'$\displaystyle \chi \approx %4.2g$'%results['chi'],bbox={'facecolor':'white','alpha':0.9,'pad':5})
                ax3.set_xlim([0,rho_max])
                ax4.set_xlim([0,rho_max])
            
                
            ax2.set_xlim([-0.5,0.1])
            ax1.set_xlim([-0.5,0.1])
            #ax2.set_xlim([-0.5,0.2])
            ax1.set_ylim([0,y_border])
                
    #mpl.rcParams['xtick.labelsize'] = 20
    #ax2.rcParams['xtick.labelsize'] = 10
    #plt.rcParams['ytick.labelsize'] = 20
    ##ax2.rcParams['ytick.labelsize'] = 20  
    #plt.rcParams['font.size'] = 12
    
    if save:
        sv_name = './../paper draft/inhib only/pics/distr.%s' % file_format
        plt.savefig(sv_name,format=file_format,dpi=600)
        print 'Figure saved as "%s"' % sv_name
        
    #ax4.plot(results['I_range'],results['f_I'])
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95, wspace=0.4, hspace=0.4)
    plt.show(block=False)