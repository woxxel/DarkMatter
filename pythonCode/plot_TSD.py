import imp

imp.load_source('sharks', 'sharkfins.py')
from sharks import *
imp.load_source('single', 'get_samples.py')
from single import *


def plot_TSD(steps=1000,rateWnt=[[1]],alpha_0=[[0,0.02,0.04]],tau_G=[[30]],eps=[[0]],eta=[[0.9]],n=[[0]],Npop=2,drive=0,save=0,file_format='png'):

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10

    fig, ax = plt.subplots(1,3,figsize=(7.5,2))

    vborder = 0.21
    hborder = 0.1
    vspace = 0.1
    hspace = 0.07

    vsize = 0.6
    hsize = 0.24

    #ax[1].set_position([0.415,0.21,0.22,0.6])
    #ax[2].set_position([0.76,0.21,0.22,0.6])

    big_ax = plt.axes([0.05,0.12,0.8,0.8])
    # big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['left'].set_visible(False)
    plt.setp(big_ax,xticks=[],yticks=[])

    big_ax.set_xlabel(r'$\displaystyle \bar{\nu}\,$[Hz]')
    big_ax.set_ylabel(r'$\displaystyle I$')

    title_descr = 0
    title_x = 0
    x_ax_string = 'rateWnt'

    axLab = ['a)','b)','c)']
    n = [[[0.0]],[[0.2]],[[0.5]]]

    for i in range(len(n)):
        results_stats = sharkfins(steps=steps,rateWnt=rateWnt,alpha_0=alpha_0,tau_G=tau_G,n=n[i],eps=eps,eta=eta,tau_M=10,Npop=Npop,drive=0,mode_calc='exact',mode_stats=1)

        ax[i].set_position([hborder+i*(hsize+hspace),vborder,hsize,vsize])

        ax[i].spines['right'].set_color('none')
        ax[i].yaxis.set_ticks_position('left')
        ax[i].spines['top'].set_color('none')
        ax[i].xaxis.set_ticks_position('bottom')

        nu_imp_idx = np.zeros(len(results_stats['gamma']))
        try:
            nu_imp_idx[0] = np.where(results_stats[x_ax_string][0]==results_stats['nu_implausible'][0])[0][0]
        except:
            nu_imp_idx[0] = -1

        #print results_stats['alpha'][
        #print results_stats[x_ax_string]

        ax[i].plot(results_stats[x_ax_string][0],results_stats['alpha'][0],'r:')
        if (len(results_stats[x_ax_string]) > 1):
            ax[i].plot(results_stats[x_ax_string][2],results_stats['alpha'][2],'r:')#,label=r'$\displaystyle \alpha$')
        ax[i].plot(results_stats[x_ax_string][0][nu_imp_idx[0]:],results_stats['sigma_V'][0][nu_imp_idx[0]:],'r:')

        ax[i].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha'][0][:nu_imp_idx[0]],'r-')#,label='inh.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
        ax[i].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['alpha_exc'][0][:nu_imp_idx[0]],'k-')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')

        #if (len(results_stats[x_ax_string]) > 1):
            #ax[0,1].plot(results_stats[x_ax_string][2][:nu_imp_idx[2]],results_stats['alpha'][2][:nu_imp_idx[2]],'r-')#,label=r'$\displaystyle \alpha$')
        ax[i].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V'][0][:nu_imp_idx[0]],'r--')#,label=r'$\displaystyle \sigma_V$')
        ax[i].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['sigma_V_exc'][0][:nu_imp_idx[0]],'k--')#,label=r'$\displaystyle \sigma_V$')

        #ax[0,1].plot(results_stats[x_ax_string][0][:nu_imp_idx[0]],results_stats['I_balance'][0][:nu_imp_idx[0]],'-',color=[0.7,0.7,0.7],label=r'$\displaystyle I_{balance}$')

        #if (len(results_stats[x_ax_string]) > 1):
            #if not np.isnan(results_stats['nu_c'][2]):
                #ax[0,2].plot([results_stats['nu_c'][2],results_stats['nu_c'][2]],[0,1],'k:')


        #ax[i].text(0.2,0.01,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)
        #if max(n)>0:
            #ax[i].text(0.05,0.07,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k,A}^2 + \sigma_{V_k,N}^2 + \sigma_{V_k,G}^2}$',fontsize=10)
        #else:
            #ax[i].text(0.075,0.043,r'$\displaystyle \sigma_{V_k} = \sqrt{\sigma_{V_k,A}^2 + \sigma_{V_k,G}^2}$',fontsize=10)

        #ax[0,1].text(0.05,0.016,r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$',fontsize=10)
        if (len(results_stats[x_ax_string]) > 1):
            if not np.isnan(results_stats['nu_c'][2]):
                ax[i].text(results_stats['nu_c'][2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)


        #ax[0,1].set_ylim([0.0,0.05])
        #plt.setp(ax[0,1],yticks=np.linspace(0.02,0.05,4),xticks=np.linspace(0,0.8,5))
        #ax[i].set_ylabel(r'$\displaystyle I$')

        #if x_ax_string == 'rateWnt':
        plt.setp(ax[i],xticks=np.linspace(0,10,6),yticks=np.linspace(0,0.15,4),xlim=[0,results_stats[x_ax_string][0][-1]],ylim=[0,0.15])
        #else:
            #plt.setp(ax[0,i],xticks=np.linspace(0,1,6),yticks=np.linspace(0,0.1,3),xlim=[0,results_stats[x_ax_string][0][-1]])
        #print "lim: ", results_stats[x_ax_string][0][-1]
        #ax[0,1].set_xlim([0,results_stats[x_ax_string][0][-1]])
        #ax[0,1].legend(prop={'size':10},loc='lower right')
        #ax[0,2].set_yticks(np.linspace(0,0.15,4))
        #ax[0,2].set_yticklabels(np.linspace(0,0.15,4),fontsize=12)
        if title_descr:
            ax[i].set_title('b) variances',position=(0.1,1.05))#,loc='left')
        else:
            ax[i].set_title('$\displaystyle b = %3.1g$' %n[i][0][0],fontsize=12)#,position=(title_x,1.05))

    ax[2].plot([0,0],[0,0],'k-',label='exc.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')
    ax[2].plot([0,0],[0,0],'r-',label='inh.')#,label=r'$\displaystyle \alpha = \sqrt{\alpha_I^2 + \alpha_0^2}$')

    ax[0].plot([0,0],[0,0],'k--',label=r'$\displaystyle \sigma_V$')
    ax[0].plot([0,0],[0,0],'k-',label=r'$\displaystyle \alpha$')

    ax[0].legend(prop={'size':10},loc='lower right')
    ax[2].legend(prop={'size':10},loc='upper left')

    if save:
        #if Npop == 1:
            #sv_name = './../paper draft/inhib only/pics/threshold_stats_%s_drive_%d.%s' % (x_ax_string,drive,file_format)
        #if Npop == 2:
        sv_name = './../paper draft/two_pop/pics/variance_stats.%s' % (file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)
