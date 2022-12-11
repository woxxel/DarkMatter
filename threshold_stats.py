import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

def threshold_stats(steps=1000,rateWnt=[0,20],alpha_0=[0.02,0.04],tau_G=[0.005],eps=[0.5],eta=[0.9],n=[0],J=-1.,Npop=1,drive=0,save=0,file_format='png',rerun=False,compile=False):

    options = {
        'Npop': Npop,
        'order': ['rateWnt','alpha_0','tau_G','n','eta','eps'],
        'rateWnt': rateWnt,
        'alpha_0': alpha_0,
        'tau_G': tau_G,
        'eps': eps,
        'eta': eta,
        'n': n,
        'drive': drive,
        'mode_stats': 1,
        'J': J
    }

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    # return res

    if len(rateWnt) > 1:
        x_key = 'rateWnt'
    if len(eps) > 1:
        x_key = 'eps'
    if len(eta) > 1:
        x_key = 'eta'
    if len(n) > 1:
        x_key = 'n'
    if len(tau_G) > 1:
        x_key = 'tau_G'

    currents = True

    set_plot_params()

    #mpl.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?

    plt_para = {
        'title': {
            'descr': False,
            'x_offset': -0.2
            },
        'two_row': True
    }

    if plt_para['two_row']:
        fig, ax = plt.subplots(2,3,figsize=(7.5,4),dpi=300)
    else:
        fig, ax = plt.subplots(2,3,figsize=(7.5,2))
        fig.delaxes(ax[1,0])
        fig.delaxes(ax[1,1])
        fig.delaxes(ax[1,2])
        ax[0,0].set_position([0.1,0.25,0.22,0.6])
        ax[0,1].set_position([0.415,0.25,0.22,0.6])
        ax[0,2].set_position([0.76,0.25,0.22,0.6])



    #big_ax = fig.add_subplot(111)
    big_ax = plt.axes([0.1,0.05,0.8,0.8])
    big_ax.set_facecolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])

    if x_key == 'rateWnt':
        big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    elif x_key == 'eps':
        big_ax.set_xlabel(r'$\displaystyle \varepsilon$')
    elif x_key == 'n':
        big_ax.set_xlabel(r'$\displaystyle b$')

    #print res['trans_implausible']
    #print res
    steps1 = res['gamma'].shape[1]
    pCol = ['r','k']

    ## define dictionary with transition point indices
    trans_idx = {}
    for key in ['inc','imp','DM','np']:
        trans_idx[key] = np.zeros(steps1,'int')
        for a in range(steps1):
            idx = np.where(res[x_key]==res[key+'_trans'][0,a])[0]
            trans_idx[key][a] = idx[0] if len(idx) else -1

    x_lim = res['inc_trans'][0,0]
    plot_q(ax[0,0],res,x_key,trans_idx,plt_para,x_lim,order=0)
    # ax[0,0].text(8,43,r'$\displaystyle q\approx\bar{\nu}^2$',fontsize=10)
    ax[0,0].legend(prop={'size':10},bbox_to_anchor=(0.5,0.05),loc='lower left',handlelength=1)

    plot_q_zoom(ax[0,1],res,x_key,trans_idx,plt_para,1,order=0)
    # ax[0,1].text(0.1,3.5,r'$\displaystyle q\approx \frac{\bar{\nu}\nu_{max}}{\sqrt{2}}$',fontsize=10)
    plot_currents(ax[0,2],res,x_key,trans_idx,plt_para,x_lim,order=0)
    plot_gamma(ax[1,0],res,x_key,trans_idx,plt_para,x_lim,order=0)
    plot_chi(ax[1,1],res,x_key,trans_idx,plt_para,x_lim,order=0)

    set_title(ax[0,0],1,'',(-0.075,0),10)
    set_title(ax[0,1],2,'',(-0.075,0),10)
    set_title(ax[0,2],3,'',(-0.075,0),10)
    set_title(ax[1,0],4,'',(-0.075,0),10)
    set_title(ax[1,1],5,'',(-0.075,0),10)
    set_title(ax[1,2],6,'',(-0.075,0),10)

    plt.setp(ax[0,0],xlabel='')
    plt.setp(ax[0,1],xlabel='')
    plt.setp(ax[0,2],xlabel='')
    plt.setp(ax[1,0],xlabel='')
    plt.setp(ax[1,1],xlabel='')

    options = {
        'rateWnt': [0,20],
        'alpha_0': [0,0.2],
        'tau_G': [0.005],
        'J': J
    }
    results_bounds = darkMatter(steps=500,options=options,rerun=rerun)
    plot_regions(ax[1,2],results_bounds,trans_idx,plt_para,res[x_key][-1],order=0)
    plt.setp(ax[1,2],xlabel='')

    plt.subplots_adjust(left=0.075, bottom=0.11, right=0.99, top=0.925, wspace=0.4, hspace=0.6)
    for i in range(2):
        for j in range(3):
            ax[i,j].spines['right'].set_color('none')
            ax[i,j].yaxis.set_ticks_position('left')
            ax[i,j].spines['top'].set_color('none')
            ax[i,j].xaxis.set_ticks_position('bottom')

    # fig.tight_layout(pad=0.5)
    if save:
        sv_name = './figures/heterogeneity.%s' % (file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)

    return res

    #if x_key == 'rateWnt':
    if currents:
    #elif x_key in ['eps','eta','n','tau_G']:
        #elif x_key in ['eps','eta','n']:
        #ax[0,0].plot(res[x_key][0][:trans_idx['imp'][0]],res['q'][0][:trans_idx['imp'][0]],'r',label='inh.')
        #ax[0,0].plot(res[x_key][0][:trans_idx['imp'][0]],res['q_exc'][0][:trans_idx['imp'][0]],'k',label='exc.')#r'$\displaystyle \alpha_0 = 0$')

        #ax[0].legend(prop={'size':10},loc='lower left')

        ax[0,0].legend(prop={'size':10},loc='lower left')
        #ax[0,0].set_ylim([2,4])
        #plt.setp(ax[0,0],yticks=np.linspace(2,4,5))

        tau_M = 0.010

        tau_A = 0.005
        #tau_G = 0.030
        tau_N = 0.200

        J_EE = res['eta'] * res['eps']*tau_M
        J_EI = np.sqrt(1-(res['eta']*res['eps'])**2)*tau_M


        var_V_A = J_EE**2 * res['rateWnt'][0] / (tau_A + tau_M) * ( (1-res['n'][0])**2/2 + (1-res['n'][0])*res['n'][0]*tau_A / (tau_A + tau_N) )
        #from excitatory NMDA synapses
        var_V_N = J_EE**2 * res['rateWnt'][0] / (tau_N + tau_M) * (res['n'][0]**2/2 + (1-res['n'][0])*res['n'][0]*tau_N / (tau_A + tau_N) );
        #from inhibitory GABA synapses
        var_V_G = J_EI**2 * res['rateWnt'][0] * 0.5 / (res['tau_G'][0] + tau_M)

        #return res, var_V_A, var_V_G
        #print res['n']
        #print var_V_N
        #print var_V_G
        ax[0,0].plot(res[x_key][0],np.sqrt(var_V_G[0]),'k--',label=r'$\displaystyle \sigma_{V_E^G}$ (GABA)')
        ax[0,0].plot(res[x_key][0],np.sqrt(var_V_A[0]),'k:',label=r'$\displaystyle \sigma_{V_E^A}$ (AMPA)')
        #ax[0,0].plot(res[x_key][0],np.sqrt(var_V_G[0]+var_V_A[0]),'k-')
        if max(n) > 0:
            ax[0,0].plot(res[x_key][0],np.sqrt(var_V_N[0]),'k-.',label=r'$\displaystyle \sigma_{V_E^N}$ (NMDA)')
            ax[0,0].legend(prop={'size':10},bbox_to_anchor=(0.05,1.2),loc='upper left',ncol=1)
        else:
            ax[0,0].text(0.2,0.005,r'$\displaystyle \sigma_{V_{EE}}$ (AMPA) $\displaystyle \hat{=} \sigma_{V_E^A}$',fontsize=10)
            ax[0,0].text(0.1,0.038,r'$\displaystyle \sigma_{V_{EI}}$ (GABA) $\displaystyle \hat{=} \sigma_{V_E^G}$',fontsize=10)

        # print(res[x_key])
        # print(res['q'])
        if x_key == 'rate':
            plt.setp(ax[0,0],xticks=np.linspace(0,10,6),yticks=np.linspace(0,0.15,4),xlim=[0,res[x_key][-1]],ylim=[0,0.15])
        else:
            #plt.setp(ax[0,0],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,4),xlim=[0,res[x_key][0][-1]])#,ylim=[0,0.075])
            plt.setp(ax[0,0],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,6),xlim=[0,res[x_key][-1]])
        ax[0,0].set_ylabel(r'$\displaystyle I$')
        ax[0,0].set_title(r'a)',position=(title_x,1.05))





    #if x_key == 'rateWnt':
    if currents:
    #elif x_key in ['eps','eta','n','tau_G']:
        # if (len(res[x_key]) > 1):
        #     ax[0,1].plot(res[x_key][2],res['alpha'][0,2],'r:')#,label=r'$\displaystyle \alpha$')
        for p in range(Npop):
            ax[0,1].plot(res[x_key][trans_idx['imp'][0]:],res['alpha'][p,0,trans_idx['imp'][0]:],pCol[p]+':')
            ax[0,1].plot(res[x_key][:trans_idx['imp'][0]],res['alpha'][p,0,:trans_idx['imp'][0]],pCol[p]+'-',label='inh.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
            ax[0,1].plot(res[x_key][trans_idx['imp'][0]:],res['sigma_V'][p,0,trans_idx['imp'][0]:],pCol[p]+':')
            ax[0,1].plot(res[x_key][:trans_idx['imp'][0]],res['sigma_V'][p,0,:trans_idx['imp'][0]],pCol[p]+'--')#,label=r'$\displaystyle \sigma_V$')

        #if (len(res[x_key]) > 1):
            #ax[0,1].plot(res[x_key][2][:trans_idx['imp'][2]],res['alpha'][2][:trans_idx['imp'][2]],'r-')#,label=r'$\displaystyle \alpha$')
        #ax[0,1].plot(res[x_key][0][:trans_idx['imp'][0]],res['I_balance'][0][:trans_idx['imp'][0]],'-',color=[0.7,0.7,0.7],label=r'$\displaystyle I_{balance}$')

        if (steps > 1):
            if not np.isnan(res['inc_trans'][0,2]):
                ax[0,2].plot([res['inc_trans'][0,2],res['inc_trans'][0,2]],[0,1],'k:')


        ax[0,1].text(0.1,0.015,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)
        if x_key == 'n':
            ax[0,1].text(0.05,0.05,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k^A}^2 + \sigma_{V_k^N}^2 + \sigma_{V_k^G}^2}$',fontsize=10)
        else:
            ax[0,1].text(0.05,0.043,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k^A}^2 + \sigma_{V_k^G}^2}$',fontsize=10)

        #ax[0,1].text(0.05,0.016,r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$',fontsize=10)
        if (len(res[x_key]) > 1):
            if not np.isnan(res['inc_trans'][0,2]):
                ax[0,1].text(res['inc_trans'][0,2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)


        #ax[0,1].set_ylim([0.0,0.05])
        #plt.setp(ax[0,1],yticks=np.linspace(0.02,0.05,4),xticks=np.linspace(0,0.8,5))
        ax[0,1].set_ylabel(r'$\displaystyle I$')

        if x_key == 'rateWnt':
            plt.setp(ax[0,1],xticks=np.linspace(0,10,6),yticks=np.linspace(0,0.15,4),xlim=[0,res[x_key][-1]],ylim=[0,0.15])
        else:
            #plt.setp(ax[0,1],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,4),xlim=[0,res[x_key][0][-1]])
            plt.setp(ax[0,1],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.075,6),xlim=[0,res[x_key][-1]])
        #print "lim: ", res[x_key][0][-1]
        #ax[0,1].set_xlim([0,res[x_key][0][-1]])
        #ax[0,1].legend(prop={'size':10},loc='lower right')
        #ax[0,2].set_yticks(np.linspace(0,0.15,4))
        #ax[0,2].set_yticklabels(np.linspace(0,0.15,4),fontsize=12)
        if title_descr:
            ax[0,1].set_title('b) variances',position=(0.1,1.05))#,loc='left')
        else:
            ax[0,1].set_title(r'b)',position=(title_x,1.05))





    #ax[0,2].plot(res[x_key][0],np.zeros(len(res[x_key][0])),'k:',linewidth=0.5)

    #if x_key == 'rateWnt':

    pLabels = ['inh.','exc.']
    #if x_key in ['eps','eta','n','tau_G']:
    if currents:
        for p in range(Npop):
            ax[0,2].plot(res[x_key][:trans_idx['imp'][0]],-res['I_balance'][p,0,:trans_idx['imp'][0]],'-',color=pCol[p],label=pLabels[p])

        if x_key == 'rateWnt':
            plt.setp(ax[0,2],xticks=np.linspace(0,10,6),yticks=np.linspace(-0.15,0,4),xlim=[0,res[x_key][-1]],ylim=[-0.15,0])
        else:
            plt.setp(ax[0,2],xticks=np.linspace(0,1,6),yticks=np.linspace(-0.15,0,4),xlim=[0,res[x_key][-1]])
        ax[0,2].legend(prop={'size':10},bbox_to_anchor=(0.95,1.1),loc='upper right')







        #ax[i,j].tick_params(axis='both', which='major', labelsize=16)
        #ax[i,j].tick_params(axis='both', which='minor', labelsize=16)
        #ax[i,j].set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    #mpl.rcParams['xtick.labelsize'] = 12
    #plt.rcParams['xtick.labelsize'] = 'x-large'

    if save:
        if Npop == 1:
            sv_name = './../paper draft/inhib only/pics/threshold_stats_%s_drive_%d.%s' % (x_key,drive,file_format)
        if Npop == 2:
            sv_name = './../paper draft/two_pop/pics/threshold_stats_%s_drive_%d.%s' % (x_key,drive,file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)
    return res

    #print res.keys()
    #y = res['sigma_V'][0]/res['I_balance'][0]
    #print res['sigma_V'][0]
    #print res['I_balance'][0]
    #print y
    #plt.figure()
    #plt.plot(res[x_key][0],y)
    #plt.plot(res[x_key][0],res['alpha_raw'][0]/max(res['alpha'][0]),'r--')
    #plt.plot(res[x_key][0],res['alpha'][0]/max(res['alpha'][0]),'r-')

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
    print("v plots: ", v_plots)

    big_ax1 = plt.axes([0.1,0.05,0.35,0.9])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    remove_frame(big_ax1)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_xlabel(r'$\displaystyle I$')

    big_ax1 = plt.axes([0.05,0.05,0.35,0.75 ])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    remove_frame(big_ax1)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_ylabel(r'$\displaystyle \rho(I)$',color='r')

    big_ax1 = plt.axes([0.05,0.25,0.35,0.75])
    big_ax1.set_axis_bgcolor('none')
    big_ax1.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    remove_frame(big_ax1)
    plt.setp(big_ax1,xticks=[],yticks=[])
    big_ax1.set_ylabel(r'$\displaystyle \nu$ [Hz]')

    big_ax2 = plt.axes([0.6,0.075,0.35,0.9])
    big_ax2.set_axis_bgcolor('none')
    big_ax2.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    remove_frame(big_ax2)
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
                print(nu_max)
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
        print('Figure saved as "%s"' % sv_name)

    #ax4.plot(results['I_range'],results['f_I'])
    #plt.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95, wspace=0.4, hspace=0.4)
    plt.show(block=False)
