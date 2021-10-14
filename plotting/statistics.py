import matplotlib.pyplot as plt
import numpy as np


def plot_q(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',approx=False):

    Npop,steps1,steps = res['q'].shape

    ## plotting q
    ax.plot(res[x_key],res['q'][0,0,:],'k:')
    ax.plot(res[x_key][:trans_idx[bound][0]],res['q'][0,0,:trans_idx[bound][0]],'k',label=r'$\displaystyle \alpha_0 = 0$')
    if (steps1 > 1):

        ax.plot(res[x_key],res['q'][0,2,:],'r:')
        ax.plot(res[x_key][:trans_idx[bound][2]],res['q'][0,2,:trans_idx[bound][2]],'r',label=r'$\displaystyle \alpha_0 = 0.04$')

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.1*steps)],res[x_key][:int(0.1*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        print('plot approx')
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--',label=r'$\displaystyle \alpha_0 = 0$')
        if (steps1 > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--',label=r'$\displaystyle \alpha_0 = 0.04$')


    ax.text(8,43,r'$\displaystyle q\approx\bar{\nu}^2$',fontsize=10)

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,200],
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
    )

    ax.legend(prop={'size':10},bbox_to_anchor=(0,1.2),loc='lower right')
    if plt_para['title']['descr']:
        ax.set_title('a) second moment',position=(0.1,1.05))#,loc='left'
    else:
        ax.set_title(r'a)',position=(plt_para['title']['x_offset'],1.05))



def plot_q_zoom(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',approx=False):

    Npop,steps1,steps = res['q'].shape

    ## plotting q
    ax.plot(res[x_key],res['q'][0,0,:],'k')
    # if (len(res[x_key]) > 1):
    for a in range(steps1):
        col = a/float(steps1)
        ax.plot(res[x_key],res['q'][0,a,:],color=(col,0,0))

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.012*steps)],res[x_key][:int(0.012*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)

        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (len(res[x_key]) > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    ax.text(0.1,3.5,r'$\displaystyle q\approx \frac{\bar{\nu}\nu_{max}}{\sqrt{2}}$',fontsize=10)

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,5],
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
    )

    if plt_para['title']['descr']:
        ax.set_title(r'b) the low $\displaystyle \bar{\nu}$ limit',position=(0.1,1.05))#loc='left')
    else:
        ax.set_title(r'b)',position=(plt_para['title']['x_offset'],1.05))


def plot_currents(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',approx=False):

    Npop,steps1,steps = res['q'].shape

    ## plot constant lines
    ax.plot([0,x_lim],[0,0],'k--',lw=0.5)
    # for a in range(steps1):
        # col = a/float(steps1)
    if (steps1 > 1):
        if not np.isnan(res['DM_trans'][0,2]):  # only possible, when there is a DM transition
            ax.axvline(res['DM_trans'][0,2],color='k',ls=':')
            ax.text(res['DM_trans'][0,2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)

    ## plotting currents
    ax.plot(res[x_key],-res['I_balance'][0,0,:],':',color=[0.7,0.7,0.7])
    ax.plot(res[x_key][:trans_idx[bound][0]],-res['I_balance'][0,0,:trans_idx[bound][0]],'-',color=[0.7,0.7,0.7],label=r'$\displaystyle I_{balance}$')

    ## plotting temporal fluctuations
    ax.plot(res[x_key][:trans_idx[bound][0]],res['sigma_V'][0,0,:trans_idx[bound][0]],'k--',label=r'$\displaystyle \sigma_V$')
    ax.plot(res[x_key][trans_idx[bound][0]:],res['sigma_V'][0,0,trans_idx[bound][0]:],'k:')

    ## plotting quenched fluctuations
    ax.plot(res[x_key][:trans_idx[bound][0]],res['alpha'][0,0,:trans_idx[bound][0]],'k-',label=r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
    ax.plot(res[x_key][trans_idx[bound][0]:],res['alpha'][0,0,trans_idx[bound][0]:],'k:')

    if (steps1 > 1):
        ax.plot(res[x_key],res['alpha'][0,2,:],'r:')#,label=r'$\displaystyle \alpha$')
        ax.plot(res[x_key][:trans_idx[bound][2]],res['alpha'][0,2,:][:trans_idx[bound][2]],'r-')#,label=r'$\displaystyle \alpha$')

    ax.text(6,-0.23,r'$\displaystyle \bar{I}_0-\Psi_0$',fontsize=10)
    ax.text(6,0.17,r'$\displaystyle \sigma_{V_k}$',fontsize=10)
    ax.text(8,0.02,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)
    #ax.set_xlim([0,15])
    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[-0.3,0.3],
        ylabel=r'$\displaystyle \bar{I} - \Psi_0$'
    )

    if plt_para['title']['descr']:
        ax.set_title('c) variances',position=(0.1,1.05))#,loc='left')
    else:
        ax.set_title(r'c)',position=(plt_para['title']['x_offset'],1.05))


def plot_gamma(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',approx=False):

    Npop,steps1,steps = res['q'].shape

    ax.plot(res[x_key],np.ones(steps),'k:',linewidth=0.5)

    for a in range(steps1):
        col = a/float(steps1)
        # if a!=1:
        for p in range(Npop):
            ax.plot(res[x_key],res['gamma'][p,a,:]**2,color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx[bound][a]],res['gamma'][p,a,:trans_idx[bound][a]]**2,color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['gamma_approx'][p,a,:]**2,color=(col,0,0),ls=':')

    if plt_para['title']['descr']:
        ax.set_title(r'e)$\displaystyle\quad$ DM exponent $\displaystyle \gamma$',position=(0.1,1.05))#,loc='left')
    else:
        ax.set_title(r'e)',position=(plt_para['title']['x_offset'],1.05))

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,7],
        ylabel=r'$\displaystyle \gamma^2$'
    )


def plot_chi(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',approx=False):

    Npop,steps1,steps = res['q'].shape

    for a in range(steps1):
        col = a/float(steps1)
        # if a!=1:
        for p in range(Npop):
            ax.plot(res[x_key],res['chi'][p,a,:],color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx[bound][a]],res['chi'][p,a,:trans_idx[bound][a]],color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['chi_approx'][p,a,:]**2,color=(col,0,0),ls=':')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[-0.5,2],
        ylabel=r'$\displaystyle \chi$'
    )

    if plt_para['title']['descr']:
        ax.set_title(r'd) skewness coeff. $\displaystyle \chi$',position=(0.125,1.05))#,loc='left')#,fontsize=12)
    else:
        ax.set_title(r'd)',position=(plt_para['title']['x_offset'],1.05))

def plot_regions(ax,res,trans_idx,plt_para,x_lim):

    Npop,steps1,steps = res['q'].shape

    mask_DM = ~np.isnan(res['DM_trans'][0,:])
    mask_no_peak = ~np.isnan(res['np_trans'][0,:])
    mask_inc = ~np.isnan(res['inc_trans'][0,:])

    ax.plot(res['DM_trans'][0,mask_DM],res['alpha_0'][mask_DM],'r-',label=r'$\displaystyle \bar{\nu}_{DM}$')
    ax.plot(res['np_trans'][0,mask_no_peak],res['alpha_0'][mask_no_peak],'k--',label=r'$\displaystyle \bar{\nu}_{no\,peak}$')
    ax.plot(res['inc_trans'][0,mask_inc],res['alpha_0'][mask_inc],'k',label=r'$\displaystyle \bar{\nu}_{inc}$')
    ax.plot(res['imp_trans'][0,mask_inc],res['alpha_0'][mask_inc],'k:',label=r'$\displaystyle \bar{\nu}_{imp}$')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,0.2],
        ylabel=r'$\displaystyle \alpha_0$'
    )

    # ax.set_yticks(np.linspace(0,0.2,6))

    ax.legend(prop={'size':10},ncol=2,bbox_to_anchor=(-0.05,1.3),loc='upper left')
    if plt_para['title']['descr']:
        ax.set_title(r'f) boundaries',position=(0.1,1.05))#,loc='left')
    else:
        ax.set_title(r'f)',position=(plt_para['title']['x_offset'],1.05))

    # if x_key in ['eps','eta','n','tau_G']:
    #     #print "eps: ", res['eps'][0]
    #     #print "eta: ", res['eta'][0][0]
    #     I_I_per_nu = np.sqrt(1-res['eps'][0]**2) - res['eps'][0]
    #
    #     eta = [0.9,0.6,0.2]
    #     for i in range(len(eta)):
    #         I_E_per_nu = np.sqrt(1-(eta[i]*res['eps'][0])**2) - eta[i]*res['eps'][0]
    #         col = float(i)/len(eta)
    #         ax.plot(res[x_key][0],I_E_per_nu,color=(col,col,col),label=r'$\displaystyle \eta = %3.1g$'%eta[i])
    #
    #     ax.plot(res[x_key][0],I_I_per_nu,'r')
    #     ax.legend(prop={'size':10},loc='lower left')
    #     ax.set_ylabel(r'$\displaystyle I^{ext} / \bar{\nu}$')
    #
    #     plt.setp(ax,xticks=np.linspace(0,0.8,5),yticks=np.linspace(0,1,3),xlim=[0,res[x_key][0][-1]])
    #
    #     ax.set_title(r'f)',position=(title_x,1.05))
