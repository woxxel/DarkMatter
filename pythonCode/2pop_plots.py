import imp

imp.load_source('sharks', 'sharkfins.py')
from sharks import *

def plot_2pop(steps=100,rate=[[0,16]],alpha_0=[[0,0.12]],tau_G=[[30]],n=[[0]],eps=None,eta=None,tau_A=5,kappa=1,drive=0,mode_calc='exact',save=0):
    
    
    if len(n[0]) > 1:
        n_plot = True
    
    para = {}
    para['n'] = n
    para['alpha_0'] = alpha_0
    para['eps'] = eps
    para['eta'] = eta
    para['tau_G'] = tau_G
    para['rateWnt'] = rate
    
    #const_labels = []
    para_len = 1
    #const_label = ''
    #const_label2 = ''
    para_labels = ['alpha_0','n','eps','eta','tau_G','rateWnt']
    idx = 0
    for key in para_labels:
        #print para[key]
        if (len(para[key]) > 1):
            #const_labels = key
            para_len = len(para[key])
            if key == 'n':
                sv_label = 'n'
                const_label = r'$\displaystyle n = $'
                const_label2 = ''
            if key == 'alpha_0':
                sv_label = 'alpha_0'
                const_label = r'$\displaystyle \alpha_0 = $'
                const_label2 = ''
            if key == 'eps':
                sv_label = 'eps'
                const_label = r'$\displaystyle \varepsilon = $'
                const_label2 = ''
            if key == 'eta':
                sv_label = 'eta'
                const_label = r'$\displaystyle \eta = $'
                const_label2 = ''
            if key == 'tau_G':
                sv_label = 'tau_G'
                const_label = r'$\displaystyle \tau_G = $'
                const_label2 = r'$\displaystyle \,$ms'
            if key == 'rateWnt':
                sv_label = 'rate'
                const_label = r'$\displaystyle \bar{\nu} = $'
                const_label2 = r'$\displaystyle \,$Hz'
            const_idx = idx
        idx += 1
        
    for key in para.keys():
        if (len(para[key]) == 1):
            para[key]*=para_len
    
    if len(rate[0]) > 1:
        x_ax_string = 'rateWnt'
        mode = 'rateWnt'
    if len(alpha_0[0]) > 1:
        x_ax_string = 'alpha_0'
        mode = 'alpha_0'
    if len(tau_G[0]) > 1:
        x_ax_string = 'tau_G'
        mode = 'tau_G'
    if len(eps[0]) > 1:
        x_ax_string = 'eps'
        mode = 'eps'
    if len(eta[0]) > 1:
        x_ax_string = 'eta'
        mode = 'eta'
    if len(n[0]) > 1:
        x_ax_string = 'b'
        mode = 'n'
    
    #mode = 'eps'
    file_format = 'png'
    #drive = 0
    
    Npop = 2
    
    if mode == 'n':
        columns = 1
    else:
        columns = 2
    
    if not (mode == 'n'):
        # set fontsize etc
        vborder = 0.12
        hborder = 0.1
        vspace = 0.1
        hspace = 0.07
        
        vsize = 0.35
        hsize = 0.21
        
        fig, ax = plt.subplots(3,2,figsize=(7.5,4))

    else:
        vborder = 0.175
        hborder = 0.09
        vspace = 0.1
        hspace = 0.07
        
        vsize = 0.72
        hsize = 0.21
        
        fig, ax = plt.subplots(3,2,figsize=(7.5,2.5))
        fig.delaxes(ax[0,1])
        fig.delaxes(ax[1,1])
        fig.delaxes(ax[2,1])
        

    #plt.figure(figsize=(10,6))
    if para_len == 1:
        for i in [1,2]:
            for p in [0,1]:
                fig.delaxes(ax[i,p])
    #else:
        #fig, ax = plt.subplots(1,2,figsize=(3,5))
    
    #plt.axes = ...               # should have two rows, one for inhibitory, one for excitatory
    
    
    levs = range(20)
    bnw = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=len(levs)-1)
    heat = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=len(levs)-1)
    bnw_regions = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)
    
    strAdd = ["_exc",""]
    
    for i in range(len(rate)):
        
        ncid, para_list, inpara, ax_list, ax_labels, simulation = sharkfins(steps=steps,rateWnt=[para['rateWnt'][i]],alpha_0=[para['alpha_0'][i]],tau_G=[para['tau_G'][i]],n=[para['n'][i]],tau_M=10,Npop=2,eps=[para['eps'][i]],eta=[para['eta'][i]],tau_A=tau_A,kappa=kappa,mode_calc=mode_calc,drive=drive,mode_stats=0)
        
        for p in range(columns):
            #print p
            regions = ncid.variables['regions'+strAdd[p]][:]#.transpose(ax_list)[0,0,:,1:,0]
            gamma = ncid.variables['gamma'+strAdd[p]][:]
            chi = ncid.variables['chi'+strAdd[p]][:]
            
            #print gamma
            #plot_matrix = np.zeros((steps,steps))
            mask_inconsistent = (regions == 3)
            mask_no_peak = (regions == 2)
            mask_implausible = (regions == 1)
            
            mask_dark_matter = (gamma**2 < 1)
            
            plot_gamma = masked_array(gamma**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
            plot_chi = masked_array(chi,mask_inconsistent + mask_no_peak + mask_dark_matter)
            plot_regions = masked_array(regions,np.invert(mask_inconsistent + mask_no_peak))
            plot_implausible = masked_array(regions,np.invert(mask_implausible))
            
            #if not multi_plot:
                #if (Npop == 1):
                    #ax = plt.axes([0.12,0.1,0.7,0.8])
                #if (Npop == 2):
            ax[i,p].set_position([hborder+i*(hspace+hsize),vborder+p*(vspace+vsize),hsize,vsize])
            
            
            ax[i,p].tick_params(axis='both', which='major', labelsize=10)
            
            results = {}
            results['DM_trans'] = ncid.variables['DM_trans'][:]
            results['no_peak_trans'] = ncid.variables['no_peak_trans'][:]
            results['inc_trans'] = ncid.variables['inc_trans'][:]
            results['nu_implausible'] = ncid.variables['nu_implausible'][:]
                
            mask_DM = ~np.isnan(results['DM_trans'])
            mask_no_peak = ~np.isnan(results['no_peak_trans'])
            mask_inc = ~np.isnan(results['inc_trans'])
            
            x_max = max(inpara[para_list[ax_list[1]]][0])#[-1]
            y_max = max(inpara[para_list[ax_list[0]]][0])#[-1]
            bnw.set_bad('k',0.)
            
            #print simulation[para_list[ax_list[1]]]
            #print simulation[para_list[ax_list[0]]]
            #print para_list[ax_list[1]]
            #print para_list[ax_list[0]]
            pgamma = ax[i,p].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_gamma,cmap=bnw,vmin=0,vmax=2)
            
            #heat.set_bad('k',0.)
            pchi = ax[i,p].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_chi,cmap=heat,vmin=0,vmax=3)
            
            pregions = ax[i,p].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
            pimplausible = ax[i,p].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.2)
            
            idx = 0
            
            #if multi_plot:
            ax[i,p].set_xlim([simulation[para_list[ax_list[1]]][0],simulation[para_list[ax_list[1]]][-1]])
            ax[i,p].set_ylim([simulation[para_list[ax_list[0]]][0],simulation[para_list[ax_list[0]]][-1]])
            #else:
                #ax.set_xlim([0,simulation[para_list[ax_list[1]]][-1]])
                #ax.set_ylim([0,simulation[para_list[ax_list[0]]][-1]])
                
            #print [para_list[ax_list[0]],para_list[ax_list[1]]]
            
            for key in para_list:
                if key in [para_list[ax_list[0]],para_list[ax_list[1]]]:
                    ticks = 11
                    if key == 'n':
                        array = np.linspace(0,1,6)
                    elif key == 'alpha_0':
                        ###print simulation[key][-1]
                        if simulation[key][-1] > 0.1:
                            ticks = 6
                        array = np.linspace(0,0.2,ticks)
                        ###print array
                    elif key == 'tau_G':
                        if simulation[key][-1] > 50:
                            ticks = 6
                        array = np.linspace(0,100,ticks)
                    elif key == 'rateWnt':
                        if simulation[key][-1] > 10:
                            ticks = 6
                        array = np.linspace(0,20,ticks)
                    elif key == 'eps':
                        array = np.linspace(0,0.8,5)
                    elif key == 'eta':
                        array = np.linspace(0,1,6)
                    
                    if idx:
                        try:
                            x_array = array[:np.where(array>simulation[key][-1])[0][0]]
                        except:
                            x_array = array[:]
                        #if multi_plot:
                        ax[i,p].set_xticks(x_array)
                            #ax.set_xticklabels(x_array,fontsize=10)
                        #else:
                            #ax.set_xticks(x_array)
                            ##ax.set_xticklabels(x_array,fontsize=12)
                    else:
                        try:
                            y_array = array[:np.where(array>simulation[key][-1])[0][0]]
                        except:
                            y_array = array[:]
                        #if multi_plot:
                        ax[i,p].set_yticks(y_array)
                            #ax.set_yticklabels(y_array,fontsize=10)
                        #else:
                            #ax.set_yticks(y_array)
                            ##ax.set_yticklabels(y_array,fontsize=12)
                    idx += 1
            
            ##DM_bound[i,0] = 0
            
            #if (multi_plot):
        #if mode == 'eps':
            #ax[i,1].set_title(r'$\displaystyle \varepsilon = %4.2g$'%eps[i][0],fontsize=12)
        #elif mode == 'eta':
            #ax[i,p].set_title(r'$\displaystyle \eta = %4.2g$'%eta[i][0],fontsize=12)
        #elif mode == 'n':
            #ax[i,p].set_title(r'$\displaystyle n = %4.2g$'%n[i][0],fontsize=12)
            #else:
        #else:
        #print ax_labels
        print "paralen: ", para_len
        if (para_len > 1):
            #p-=1
            print '%s %4.2g %s' % (const_label,para[para_labels[const_idx]][i][0],const_label2)
            ax[i,p].set_title('%s %4.2g %s' % (const_label,para[para_labels[const_idx]][i][0],const_label2),fontsize=12)
    
    print 
    # axes for additional labels
    print "mode: ", mode
    if mode == 'n':
        big_ax = plt.axes([hborder-0.02,vborder,3*hsize+2*hspace,vsize])
    else:
        big_ax = plt.axes([hborder-0.02,vborder,3*hsize+2*hspace,2*vsize+vspace])
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    big_ax.set_xlabel(ax_labels[1],fontsize=12)
    big_ax.set_ylabel(ax_labels[0],fontsize=12)
    
    if not (mode == 'n'):
        inh_ax = plt.axes([hborder-0.01,vborder+vsize+vspace,hsize,vsize])
        inh_ax.set_axis_bgcolor('none')
        inh_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        inh_ax.spines['top'].set_visible(False)
        inh_ax.spines['right'].set_visible(False)
        inh_ax.spines['bottom'].set_visible(False)
        inh_ax.spines['left'].set_visible(False)
        #inh_ax.set_xlabel(ax_labels[1],fontsize=12)
        inh_ax.set_ylabel('inhibitory',fontsize=10)
        
        exc_ax = plt.axes([hborder-0.01,vborder,hsize,vsize])
        exc_ax.set_axis_bgcolor('none')
        exc_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
        exc_ax.spines['top'].set_visible(False)
        exc_ax.spines['right'].set_visible(False)
        exc_ax.spines['bottom'].set_visible(False)
        exc_ax.spines['left'].set_visible(False)
        #exc_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]',fontsize=12)
        exc_ax.set_ylabel('excitatory',fontsize=10)
    
    
    # plot colorbars
    #gamma_bar_len = vsize
    #axcb1 = plt.axes([0.9,vborder+0.23,0.03,vsize-0.23])
    ##axcb1 = plt.axes([0.9,vborder+0.2,0.03,0.55])
    #axcb2 = plt.axes([0.9,vborder,0.03,0.17])
    
    total_height = columns*vsize+(columns-1)*vspace
    cb1_length = total_height*0.65
    cb2_length = total_height*0.25
    
    axcb1 = plt.axes([0.9,vborder+total_height*0.35,0.03,cb1_length])
    axcb2 = plt.axes([0.9,vborder,0.03,cb2_length])  
    
    
    axcb1.tick_params(axis='both', which='major', labelsize=12)
    axcb2.tick_params(axis='both', which='major', labelsize=12)
    
    plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,7))
    plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,10),ticks=np.linspace(0,1,3))
    
    ###axcb2.set_yticks(np.linspace(1,0,3))
    ###axcb2.set_yticklabels(np.linspace(1,0,3))
    axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=12)
    axcb2.set_ylabel(r'$\displaystyle \gamma^2$',fontsize=12)
    
    
    if save:
        #if mode == 'eps':
        #strAdd = '%s=%4.2g'%(sv_label,eta[0])
        #if mode == 'eta':
            #strAdd = 'eps=%4.2g'%(eps[0])
        #if mode == 'n':
            #strAdd = 'eta=%4.2g_eps=%4.2g'%(eta[0],eps[0])
            
        sv_name = '../paper draft/two_pop/pics/phase_%s_steps=%d_drive%d.%s' % (sv_label,steps,drive,file_format)
        plt.savefig(sv_name,dpi=300)
        print 'Figure saved as "%s"' % sv_name
    plt.show(block=False)


#def plot_eps():




#def plot_n():




