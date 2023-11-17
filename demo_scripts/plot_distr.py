import numpy as np
import matplotlib.pyplot as plt

from darkMatter import darkMatter
from utils.plots import *
from utils.parameters import *

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
