### plot analytic approximation

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

class network:
    
    tau_I = 0.005   # synaptic timeconstant in s
    tau_m = 0.010   # membrane timeconstant in s
    J = -1 * tau_m          # synaptic coupling strength
    
    def I_squared_nu(self, nu, q):
        return - ( self.alpha(q)**2 + self.sigma_V(nu)**2 ) * np.log( (nu/self.rate_max())**2 * (1 + (self.alpha(q) / self.sigma_V(nu))**2) )

    def I_squared_q(self, nu, q):
        return -( self.alpha(q)**2 + 1./2 * self.sigma_V(nu)**2 ) * np.log( ( q/self.rate_max()**2 )**2 * (1 + 2*(self.alpha(q) / self.sigma_V(nu))**2) )



    def alpha(self, q):
        return np.sqrt(self.J**2 * q)

    def sigma_V(self, nu):
        return np.sqrt((self.J**2 * nu) / self.tau_q())
    
    def rate_max(self):
        return (2 * math.pi * np.sqrt(self.tau_I*self.tau_m))**(-1)
    
    def tau_q(self):
        return 2 * (self.tau_I + self.tau_m)

#def plot_selfcon(nu):
    



def plt_approx_single(rate=5,steps=100,save=0,file_format='png'):
    
    fig = plt.figure(figsize=(7.5,2.5))
    
    lower_bound = 0.18
    box_height = 0.7
    box_width = 0.22
    v_space = 0.05
    h_space = 0.1
    
    ax_approx = plt.axes([0.1,lower_bound,box_width,box_height])
    ax_q = plt.axes([0.1+box_width+h_space,lower_bound,box_width+0.03,box_height-0.3])
    ax_q_zoom = plt.axes([0.1+2*box_width+h_space+0.03-0.15,lower_bound+0.5,0.125,box_height/3])
    ax_gamma = plt.axes([0.1+2*(box_width+h_space)+0.03,lower_bound,box_width,(box_height-v_space)/2])
    ax_chi = plt.axes([0.1+2*(box_width+h_space)+0.03,lower_bound+(v_space + box_height)/2,box_width,(box_height-v_space)/2])
    #ax_approx = plt.axes([0.1,0.62,0.2,0.35])
    #ax_q = plt.axes([0.43,0.62,0.22,0.35])
    #ax_q_zoom = plt.axes([0.76,0.62,0.22,0.35])
    #ax_gamma = plt.axes([0.1,0.12,0.22,0.35])
    #ax_chi = plt.axes([0.43,0.12,0.22,0.35])
    #ax_distr = plt.axes([0.1,0.12,0.415,0.35])
    #ax_KL_single = plt.axes([0.76,0.12,0.22,0.35])
    
    #ax_entropy = plt.axes([0.76,0.12,0.22,0.15])
    #ax_KL_single = plt.axes([0.76,0.28,0.22,0.15])
    #ax_bit_rate = plt.axes([0.76,0.44,0.22,0.15])
    
    
    #ax_KL_phase = plt.axes([0.7,0.15,0.27,0.35])
    
    ### load data from sketch of approximation
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['font.size'] = 10
    
    nu = rate
    
    q = np.linspace(0,4*nu**2,1001)
    
    y_range = [-0.4,0]
    
    net = network()
    
    min_idx = np.argmin(abs(np.sqrt(net.I_squared_q(nu,q))-np.sqrt(net.I_squared_nu(nu,q))))
    
    #plt.figure(figsize=(3,2.3))
    #ax = plt.axes([0.2,0.2,0.65,0.75])
    
    ax_approx.plot([nu**2,nu**2],[-1,-np.sqrt(net.I_squared_nu(nu,nu**2))],color=[0.6,0.6,0.6],linewidth=3,linestyle='--')
    ax_approx.plot(q,-np.sqrt(net.I_squared_nu(nu,np.ones(len(q))*nu**2)),'k--',label=r"approx.: $\displaystyle q=\bar{\nu}^2$")       # plot approximation
    
    ax_approx.plot(q,-np.sqrt(net.I_squared_nu(nu,q)),'k-',label=r"solution for $\displaystyle \bar{\nu}$")
    ax_approx.plot(q,-np.sqrt(net.I_squared_q(nu,q)),'-',color=(0.5,0.5,0.5),label=r"solution for $\displaystyle q$")
    
    
    #plt.plot(q,net.I_squared_nu(nu,q),'k-')
    #plt.plot(q,net.I_squared_q(nu,q),'k--')
    #plt.setp(ax_approx,xticks=np.linspace(0,100,6), xticklabels=np.linspace(0,100,6).astype('int'))
    plt.setp(ax_approx,yticks=np.linspace(y_range[0],y_range[1],5))#, yticklabels=np.linspace(y_range[0],y_range[1],5))
    ax_approx.set_ylim(y_range)
    ax_approx.set_xlabel(r'$\displaystyle q\,$[Hz$\displaystyle^2$]')
    ax_approx.set_ylabel(r'$\displaystyle I_0-\Psi_0$')
    
    ax_approx.spines['right'].set_color('none')
    ax_approx.yaxis.set_ticks_position('left')
    ax_approx.spines['top'].set_color('none')
    ax_approx.xaxis.set_ticks_position('bottom')
    
    ##print q[min_idx]
    ##print np.sqrt(net.I_squared_q(nu,q[min_idx]))
    ax_approx.annotate(r'$\displaystyle q = \bar{\nu}^2$',xy=[nu**2,-np.sqrt(net.I_squared_nu(nu,nu**2))],xytext=[nu**2-22,-0.15],arrowprops=dict(arrowstyle="->"),fontsize=10)
    ax_approx.annotate(r'$\displaystyle (\bar{\nu}^{\star},q^{\star})$',xy=[q[min_idx],-np.sqrt(net.I_squared_nu(nu,q[min_idx]))],xytext=[nu**2+15,-0.35],arrowprops=dict(arrowstyle="->"),fontsize=10)
    
    ax_approx.legend(prop={'size':10},frameon=False,loc=[0.3,0.6])
    
    #sv_name = '../paper draft/inhib only/pics/selfcon_approx.png'
    #print 'Figure saved as "%s"' % sv_name
    #plt.savefig(sv_name,dpi=300)
    
    #plt.show(block=False)
    
    
    
    
    alpha_0=[[0.0,0.02,0.04]]
    ### load data from 3 example graphs for q, gamma and chi plots
    results = sharkfins(steps=steps,rateWnt=[[0,20]],alpha_0=alpha_0,tau_G=[[5]],n=[[0],[0],[0]],tau_M=10,mode_calc='exact',mode_stats=3,save=0,file_format='png')
    
    #col = []
    i_iter = results['rateWnt'].shape[0]
    
    for i in range(i_iter):
        c = float(i)/i_iter
        #print alpha_0[0
        ax_q.plot(results['rateWnt'][0],results['q'][i],'-',color=(c,0,0),label=r'$\displaystyle \alpha_0 = %4.2f$' % alpha_0[0][i])
        ax_q.plot(results['rateWnt'][0],results['q_approx'][i],'--',color=(c,0,0))
        
        ax_q_zoom.plot(results['rateWnt'][0],results['q'][i],'-',color=(c,0,0))
        ax_q_zoom.plot(results['rateWnt'][0],results['q_approx'][i],'--',color=(c,0,0))
        
        ax_gamma.plot(results['rateWnt'][0],results['gamma'][i]**2,'-',color=(c,0,0),label=r'$\displaystyle \alpha_0 = %4.2f$' % alpha_0[0][i])
        ax_gamma.plot(results['rateWnt'][0],results['gamma_approx'][i]**2,'--',color=(c,0,0))
        
        ### plotten ab nu_DM
        ax_chi.plot(results['rateWnt'][0][results['nu_c_idx'][i,0]:],results['chi'][i][results['nu_c_idx'][i,0]:],'-',color=(c,0,0))
        ax_chi.plot(results['rateWnt'][0][results['nu_c_idx'][i,0]:],results['chi_approx'][i][results['nu_c_idx'][i,0]:],'--',color=(c,0,0))
        #ax_KL_single.plot(results['rateWnt'][0],results['KL'][i],color=(c,0,0),label=r'$\displaystyle \alpha_0 = %4.2f$' % alpha_0[0][i])
        #ax_entropy.plot(results['rateWnt'][0],results['entropy'][i],color=(c,0,0),label=r'$\displaystyle \alpha_0 = %4.2f$' % alpha_0[0][i])
        #ax_bit_rate.plot(results['rateWnt'][0],results['KL'][i]/results['entropy'][i],color=(c,0,0),label=r'$\displaystyle \alpha_0 = %4.2f$' % alpha_0[0][i])
    
    x_max = 12
    ax_q.set_xlim([0,x_max])
    ax_q_zoom.set_xlim([0,1])
    ax_gamma.set_xlim([0,x_max])
    ax_chi.set_xlim([0,x_max])
    #ax_KL_single.set_xlim([0,x_max])
    #ax_entropy.set_xlim([0,x_max])
    #ax_bit_rate.set_xlim([0,x_max])
    #ax_q.set_xlim([0,x_max])
    
    ax_q.set_ylim([0,220])
    ax_q_zoom.set_ylim([0,4])
    ax_gamma.set_ylim([0,7])
    ax_chi.set_ylim([0,2])
    #ax_KL_single.set_ylim([0,0.5])
    #ax_entropy.set_ylim([-5,5])
    #ax_bit_rate.set_ylim([-1,1])
    ax_q.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    #ax_chi.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    ax_gamma.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    #ax_KL_single.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    #ax_KL_single.set_ylabel(r'KL entropy')
    
    
    ax_q.spines['top'].set_color('none')
    ax_q.xaxis.set_ticks_position('bottom')
    ax_q_zoom.spines['top'].set_color('none')
    ax_q_zoom.xaxis.set_ticks_position('bottom')
    ax_q_zoom.spines['left'].set_color('none')
    ax_q_zoom.yaxis.set_ticks_position('right')
    
    ax_chi.spines['right'].set_color('none')
    ax_chi.yaxis.set_ticks_position('left')
    ax_chi.spines['top'].set_color('none')
    ax_chi.xaxis.set_ticks_position('bottom')
    #ax_chi.yaxis.set_label_position('left')
    ax_gamma.spines['right'].set_color('none')
    ax_gamma.spines['top'].set_color('none')
    ax_gamma.xaxis.set_ticks_position('bottom')
    ax_gamma.yaxis.set_ticks_position('left')
    #ax_KL_single.spines['right'].set_color('none')
    #ax_KL_single.spines['top'].set_color('none')
    #ax_KL_single.xaxis.set_ticks_position('bottom')
    #ax_KL_single.yaxis.set_ticks_position('left')
    
    ax_q_zoom.yaxis.set_label_position('right')
    
                
    #ax_chi.set_xticks([])
    ax_q_zoom.set_xticks(np.linspace(0,1,3))
    ax_q_zoom.set_yticks(np.linspace(0,4,5))
    ax_chi.set_xticks([])
    ax_chi.set_yticks(np.linspace(0,2,3))
    ax_gamma.set_yticks(np.linspace(0,6,4))
    ax_q.set_ylabel(r'q [Hz$\displaystyle^2$]')
    ax_q_zoom.set_ylabel(r'q [Hz$\displaystyle^2$]')
    ax_chi.set_ylabel(r'$\displaystyle \chi$')
    ax_gamma.set_ylabel(r'$\displaystyle \gamma^2$')
    
    ax_q.legend(frameon=False,prop={'size':10},bbox_to_anchor=(-0.25,1.85),loc='upper left')
    #plt.plot(simulation['rateWnt'],results['KL'][0],'r--')
    
    ax_approx.text(-20,0.025,'a)',fontsize=14)
    ax_q.text(-2,405,'b)',fontsize=14)
    ax_chi.text(-2,2.25,'c)',fontsize=14)
    
    if save:
        sv_name = './../paper draft/inhib only/pics/approx.%s' % file_format
        plt.savefig(sv_name,format=file_format,dpi=600)
        print 'Figure saved as "%s"' % sv_name
    
    plt.show(block=False)

    
    ### load whole phase plot of KL-div values
    
def plt_approx_phase(steps,save=0,file_format='png'):
    
    
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #mpl.rcParams['axes.labelsize'] = 'large'
    plt.rcParams['xtick.labelsize'] = 10
    #plt.rcParams['xtick.major.size'] = 2
    plt.rcParams['ytick.labelsize'] = 10
    #plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    
  
  
    results1,simulation = sharkfins(steps=steps,rateWnt=[[0,20]],alpha_0=[[0,0.2]],tau_G=[[5]],n=[[0]],tau_M=10,mode_calc='exact',mode_stats=4,save=0,file_format='png')
    
    ### set colormaps
    levs = range(20)
    bnw = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=len(levs)-1)
    bnw.set_bad('k',0.)

    heat = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=len(levs)-1)
    heat.set_bad('k',0.)

    bnw_regions = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)
    
    cb_err = mcolors.LinearSegmentedColormap.from_list(name='err',colors=[(0,(1,1,1)),(1,(1,0,0))],N=len(levs)-1)
    
    KL_min = -3
    
    box_size = 0.22
    
    v_space = 0.09
    big_ax_v_offset = 0.025
    
    fig, ax = plt.subplots(3,3,figsize=(7.5,7.5))
    
    att_string = ['','_approx']
    
    ax[0,0].set_position([0.085,0.075+2*(box_size+v_space),box_size,box_size])
    ax[0,1].set_position([0.085+box_size+0.05,0.075+2*(box_size+v_space),box_size,box_size])
    ax[0,2].set_position([0.085+2*box_size+0.1,0.075+2*(box_size+v_space),box_size,box_size])
    
    big_ax = plt.axes([0.04,0.075+2*(box_size+v_space)-big_ax_v_offset,0.15+3*box_size,0.27])
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])
    big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    big_ax.set_ylabel(r'$\displaystyle \alpha_0$')
    
    ### prepare data
    for i in range(2):
            string = att_string[i]
            
            mask_inconsistent = (results1['regions'] == 3)
            mask_no_peak = (results1['regions'] == 2)
            mask_implausible = (results1['regions'] == 1)
            mask_dark_matter = (results1['gamma'+string]**2 < 1)
            
            plot_gamma = masked_array(results1['gamma'+string]**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
            plot_chi = masked_array(results1['chi'+string],mask_inconsistent + mask_no_peak + mask_dark_matter)
            plot_regions = masked_array(results1['regions'],np.invert(mask_inconsistent + mask_no_peak))
            plot_implausible = masked_array(results1['regions'],np.invert(mask_implausible))
            
    
            pgamma = ax[0,i].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_gamma,cmap=bnw,vmin=0,vmax=2)
            pchi = ax[0,i].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_chi,cmap=heat,vmin=0,vmax=3)
            ###bnw.set_bad((0.9,0.9,0.9),1.)
            pregions = ax[0,i].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
            pimplausible = ax[0,i].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
            
            ax[0,i].set_xlim([0,20])
            ax[0,i].set_ylim([0,0.2])
            
    
    pKL = ax[0,2].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],np.log10(results1['KL']),cmap=cb_err,vmin=KL_min,vmax=0)
    #pKL = ax[0,2].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],results1['KL']/results['entropy'],cmap=cb_err,vmin=KL_min,vmax=0)
    pregions = ax[0,2].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
    #pimplausible = ax[0,2].pcolormesh(simulation['rateWnt'],simulation['alpha_0'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
            
    ax[0,2].set_xlim([0,20])
    ax[0,2].set_ylim([0,0.2])
    
    ax[0,1].set_yticks([])
    ax[0,2].set_yticks([])
    
    ax[0,0].tick_params(axis='both', which='major', labelsize=10)
    ax[0,1].tick_params(axis='both', which='major', labelsize=10)
    ax[0,2].tick_params(axis='both', which='major', labelsize=10)
    #ax[0,2].set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    
    
    results2,simulation = sharkfins(steps=steps,rateWnt=[[0,20]],alpha_0=[[0.02]],tau_G=[[0,100]],n=[[0]],tau_M=10,mode_calc='exact',mode_stats=4,save=0,file_format='png')
    
    ax[1,0].set_position([0.085,0.075+box_size+v_space,box_size,box_size])
    ax[1,1].set_position([0.085+box_size+0.05,0.075+box_size+v_space,box_size,box_size])
    ax[1,2].set_position([0.085+2*box_size+0.1,0.075+box_size+v_space,box_size,box_size])
    
    big_ax = plt.axes([0.04,0.075+box_size+v_space-big_ax_v_offset,0.15+3*box_size,0.27])
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])
    big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    big_ax.set_ylabel(r'$\displaystyle \tau_G\,$[ms]')
    
    
    for i in range(2):
            string = att_string[i]
            
            mask_inconsistent = (results2['regions'] == 3)
            mask_no_peak = (results2['regions'] == 2)
            mask_implausible = (results2['regions'] == 1)
            mask_dark_matter = (results2['gamma'+string]**2 < 1)
            
            plot_gamma = masked_array(results2['gamma'+string]**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
            plot_chi = masked_array(results2['chi'+string],mask_inconsistent + mask_no_peak + mask_dark_matter)
            plot_regions = masked_array(results2['regions'],np.invert(mask_inconsistent + mask_no_peak))
            plot_implausible = masked_array(results2['regions'],np.invert(mask_implausible))
            
    
            pgamma = ax[1,i].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_gamma,cmap=bnw,vmin=0,vmax=2)
            pchi = ax[1,i].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_chi,cmap=heat,vmin=0,vmax=3)
            ###bnw.set_bad((0.9,0.9,0.9),1.)
            pregions = ax[1,i].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
            pimplausible = ax[1,i].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
            
            ax[1,i].set_xlim([0,20])
            ax[1,i].set_ylim([0,100])
            
    
    pKL = ax[1,2].pcolormesh(simulation['rateWnt'],simulation['tau_G'],np.log10(results2['KL']),cmap=cb_err,vmin=KL_min,vmax=0)
    pregions = ax[1,2].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
    #pimplausible = ax[1,2].pcolormesh(simulation['rateWnt'],simulation['tau_G'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
            
    ax[1,2].set_xlim([0,20])
    ax[1,2].set_ylim([0,100])
    
    ax[1,1].set_yticks([])
    ax[1,2].set_yticks([])
    
    ax[1,0].tick_params(axis='both', which='major', labelsize=10)
    ax[1,1].tick_params(axis='both', which='major', labelsize=10)
    ax[1,2].tick_params(axis='both', which='major', labelsize=10)
    
    #ax[1,2].set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    
    
    
    results3,simulation = sharkfins(steps=steps,rateWnt=[[3]],alpha_0=[[0,0.1]],tau_G=[[0,100]],n=[[0]],tau_M=10,mode_calc='exact',mode_stats=4,save=0,file_format='png')
    
    ax[2,0].set_position([0.085,0.075,box_size,box_size])
    ax[2,1].set_position([0.085+box_size+0.05,0.075,box_size,box_size])
    ax[2,2].set_position([0.085+2*box_size+0.1,0.075,box_size,box_size])
    
    big_ax = plt.axes([0.04,0.075-big_ax_v_offset,0.15+3*box_size,0.27])
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])
    big_ax.set_xlabel(r'$\displaystyle \tau_G\,$[ms]')
    big_ax.set_ylabel(r'$\displaystyle \alpha_0$')
    
    for i in range(2):
            string = att_string[i]
            
            mask_inconsistent = (results3['regions'] == 3)
            mask_no_peak = (results3['regions'] == 2)
            mask_implausible = (results3['regions'] == 1)
            mask_dark_matter = (results3['gamma'+string]**2 < 1)
            
            plot_gamma = masked_array(results3['gamma'+string]**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
            plot_chi = masked_array(results3['chi'+string],mask_inconsistent + mask_no_peak + mask_dark_matter)
            plot_regions = masked_array(results3['regions'],np.invert(mask_inconsistent + mask_no_peak))
            plot_implausible = masked_array(results3['regions'],np.invert(mask_implausible))
            
    
            pgamma = ax[2,i].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_gamma,cmap=bnw,vmin=0,vmax=2)
            pchi = ax[2,i].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_chi,cmap=heat,vmin=0,vmax=3)
            ###bnw.set_bad((0.9,0.9,0.9),1.)
            pregions = ax[2,i].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
            pimplausible = ax[2,i].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
            
            ax[2,i].set_xlim([0,100])
            ax[2,i].set_ylim([0,0.1])
            
    
    pKL = ax[2,2].pcolormesh(simulation['tau_G'],simulation['alpha_0'],np.log10(results3['KL']),cmap=cb_err,vmin=KL_min,vmax=0)
    pregions = ax[2,2].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_regions,cmap=bnw_regions,vmin=2,vmax=3)
    #pimplausible = ax[2,2].pcolormesh(simulation['tau_G'],simulation['alpha_0'],plot_implausible,cmap=bnw_regions,vmin=1,vmax=3,alpha=0.4)
    
    ax[2,2].set_xlim([0,100])
    ax[2,2].set_ylim([0,0.1])
    
    ax[2,1].set_yticks([])
    ax[2,2].set_yticks([])
    
    ax[2,0].tick_params(axis='both', which='major', labelsize=10)
    ax[2,1].tick_params(axis='both', which='major', labelsize=10)
    ax[2,2].tick_params(axis='both', which='major', labelsize=10)
    
    #ax[2,2].set_xlabel(r'$\displaystyle \tau_G\,$[ms]')
    

    axcb1 = plt.axes([0.9,0.3,0.03,0.4])
    axcb2 = plt.axes([0.9,0.075,0.03,0.2])
    axcb3 = plt.axes([0.9,0.75,0.03,0.075+3*box_size+2*v_space-0.75])
        
    axcb1.tick_params(axis='both', which='major', labelsize=12)
    axcb2.tick_params(axis='both', which='major', labelsize=12)
    axcb3.tick_params(axis='both', which='major', labelsize=12)
    
    plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,7))
    plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,10),ticks=np.linspace(0,1,3))
    plt.colorbar(pKL, cax = axcb3,boundaries=np.linspace(KL_min,0,100),ticks=np.linspace(KL_min,0,4))
    
    axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=12)
    axcb2.set_ylabel(r'$\displaystyle \gamma^2$',fontsize=12)
    axcb3.set_ylabel(r'$\displaystyle \log_{10}(KL)$',fontsize=12)
    
    #ax[0,0].text(0.22,'exact',fontsize=12)
    #ax[0,1].set_title('approximate',fontsize=12)
    #ax[0,2].set_title('KL divergence',fontsize=12)
    
    ax[0,0].set_title('exact',position=(0.5,1.2),fontsize=12)
    ax[0,1].set_title('approximate',position=(0.5,1.2),fontsize=12)
    ax[0,2].set_title('KL divergence',position=(0.5,1.2),fontsize=12)
    
    ax[0,0].text(-5,0.215,r'a) $\displaystyle \tau_G = 5\,$ms',fontsize=12)#bbox={'facecolor':'white','pad':5},fontsize=12)
    ax[1,0].text(-5,107.5,r'b) $\displaystyle \alpha_0 = 0.02$',fontsize=12)#bbox={'facecolor':'white','pad':5},fontsize=12)
    ax[2,0].text(-25,0.1075,r'c) $\displaystyle \bar{\nu} = 3\,$Hz',fontsize=12)#,bbox={'facecolor':'white','pad':5},fontsize=12)
    
    if save:
        sv_name = './../paper draft/inhib only/pics/KL_phase.%s' % file_format
        plt.savefig(sv_name,format=file_format,dpi=600)
        print 'Figure saved as "%s"' % sv_name
        
    
    plt.show(block=False)