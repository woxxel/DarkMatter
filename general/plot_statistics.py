import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.colors as mcolors

import numpy as np


def plot_q(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',idxs=None,approx=False,order=1):

    Npop,steps1,steps = res['q'].shape

    ## plotting q
    # ax.plot(res[x_key],res['q'][0,0,:],'k:')
    # ax.plot(res[x_key][:trans_idx[bound][0]],res['q'][0,0,:trans_idx[bound][0]],'k',label=r'$\displaystyle \alpha_0 = 0$')
    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        ax.plot(res[x_key],res['q'][0,i,:],color=(col,0,0),ls=':')
        ax.plot(res[x_key][:trans_idx[bound][i]],res['q'][0,i,:trans_idx[bound][i]],color=(col,0,0),label=r'$\displaystyle \alpha_0 = %g$'%res['alpha_0'][a])

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.1*steps)],res[x_key][:int(0.1*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        print('plot approx')
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (steps1 > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,200],
        # ylim=[0,20],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
    )

    if order:
        set_title(ax,order=order,title='second moment')



def plot_q_zoom(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',idxs=None,approx=False,order=1):

    Npop,steps1,steps = res['q'].shape

    ## plotting q
    ax.plot(res[x_key],res['q'][0,0,:],'k')
    # if (len(res[x_key]) > 1):
    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        ax.plot(res[x_key],res['q'][0,a,:],color=(col,0,0))

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.012*steps)],res[x_key][:int(0.012*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)

        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (len(res[x_key]) > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,5],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
    )

    if order:
        set_title(ax,order=order,title=r'the low $\displaystyle \bar{\nu}$ limit')


def plot_currents(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',idxs=None,approx=False,options=['var','I'],order=1):

    Npop,steps1,steps = res['q'].shape

    ## plot constant lines
    ax.plot([0,x_lim],[0,0],'k--',lw=0.5)

    if not idxs:
        idxs = range(steps1)

    if 'var' in options:
        if (steps1 > 1):
            if not np.isnan(res['DM_trans'][0,2]):  # only possible, when there is a DM transition
                ax.axvline(res['DM_trans'][0,2],color='k',ls=':')
                ax.text(res['DM_trans'][0,2]+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)

        ## plotting temporal fluctuations
        ax.plot(res[x_key][:trans_idx[bound][0]],res['sigma_V'][0,0,:trans_idx[bound][0]],color='k',ls='--',label=r'$\displaystyle \sigma_V$')
        ax.plot(res[x_key][trans_idx[bound][0]:],res['sigma_V'][0,0,trans_idx[bound][0]:],color='k',ls=':')

        ## plotting quenched fluctuations
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            ax.plot(res[x_key][:trans_idx[bound][i]],res['alpha'][0,i,:trans_idx[bound][i]],color=(col,0,0),ls='-',label=r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
            ax.plot(res[x_key][trans_idx[bound][i]:],res['alpha'][0,i,trans_idx[bound][i]:],color=(col,0,0),ls=':')

        ax.text(6,0.2,r'$\displaystyle \sigma_{V_k}$',fontsize=10)
        ax.text(8,0.04,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)

        plt.setp(ax,
            ylim=[-0.02,0.15],
        )
    if 'I' in options:
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            ## plotting currents
            ax.plot(res[x_key],-res['I_balance'][0,i,:],':',color=[0.7,0,0],lw=0.8)
            ax.plot(res[x_key][:trans_idx[bound][i]],-res['I_balance'][0,i,:trans_idx[bound][i]],'-',color=(col,0,0),label=r'$\displaystyle I_{balance}$',lw=0.8)

        ax.text(6,-0.22,r'$\displaystyle \bar{I}_0-\Psi_0$',fontsize=10)

        plt.setp(ax,
            ylim=[-0.3,0.02]
        )

    if 'var' in options and 'I' in options:
        plt.setp(ax,
            ylim=[-0.3,0.3]
        )


    # if (steps1 > 1):
    #     ax.plot(res[x_key],res['alpha'][0,2,:],'r:')#,label=r'$\displaystyle \alpha$')
    #     ax.plot(res[x_key][:trans_idx[bound][2]],res['alpha'][0,2,:][:trans_idx[bound][2]],'r-')#,label=r'$\displaystyle \alpha$')

    #ax.set_xlim([0,15])
    plt.setp(ax,
        xlim=[0,x_lim],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle I$'
        # ylabel=r'$\displaystyle \bar{I} - \Psi_0$'
    )

    if order:
        set_title(ax,order=order,title='variances')


def plot_gamma(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',idxs=None,approx=False,order=1):

    Npop,steps1,steps = res['q'].shape

    ax.plot(res[x_key],np.ones(steps),'k:',linewidth=0.5)

    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        # if a!=1:
        for p in range(Npop):
            ax.plot(res[x_key],res['gamma'][p,a,:]**2,color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx[bound][a]],res['gamma'][p,a,:trans_idx[bound][a]]**2,color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['gamma_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    if order:
        set_title(ax,order=order,title='$\displaystyle\quad$ DM exponent $\displaystyle \gamma$')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,7],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle \gamma^2$'
    )


def plot_chi(ax,res,x_key,trans_idx,plt_para,x_lim,bound='imp',idxs=None,approx=False,order=1):

    Npop,steps1,steps = res['q'].shape

    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        # if a!=1:
        for p in range(Npop):
            ax.plot(res[x_key],res['chi'][p,a,:],color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx[bound][a]],res['chi'][p,a,:trans_idx[bound][a]],color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['chi_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    plt.setp(ax,
        xlim=[0,x_lim],
        ylim=[0,3],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle \chi$'
    )

    if order:
        set_title(ax,order=order,title='skewness coeff. $\displaystyle \chi$')


def plot_regions(ax,res,trans_idx,plt_para,x_lim,order=1):

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
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle \alpha_0$'
    )

    # ax.set_yticks(np.linspace(0,0.2,6))

    ax.legend(prop={'size':10},ncol=2,bbox_to_anchor=(0.05,1.4),loc='upper left',handlelength=1)

    if order:
        set_title(ax,order=order,title='boundaries')

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

def plot_fins(ax,x_arr,y_arr,gamma,chi,regions,plt_para):

    levs = 20
    plt_para['cb_plotted'] = False
    plt_para['bnw'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=levs-1)
    plt_para['heat'] = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=levs-1)
    plt_para['bnw_regions'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)

    mask_inconsistent = (regions == 3)
    mask_no_peak = (regions == 2)
    mask_implausible = (regions == 1)

    mask_dark_matter = (gamma**2 < 1)

    plot_gamma = masked_array(gamma**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
    plot_chi = masked_array(chi,mask_inconsistent + mask_no_peak + mask_dark_matter)
    plot_regions = masked_array(regions,np.invert(mask_inconsistent + mask_no_peak))
    plot_implausible = masked_array(regions,np.invert(mask_implausible))

    # ax.tick_params(axis='both', which='major', labelsize=10)

    x_max = x_arr[-1]
    y_max = y_arr[-1]
    plt_para['bnw'].set_bad('k',0.)
    plt_para['heat'].set_bad('k',0.)
    pgamma = ax.pcolormesh(x_arr,y_arr,plot_gamma,cmap=plt_para['bnw'],vmin=0,vmax=2,shading='auto')
    pchi = ax.pcolormesh(x_arr,y_arr,plot_chi,cmap=plt_para['heat'],vmin=0,vmax=3,shading='auto')
    pregions = ax.pcolormesh(x_arr,y_arr,plot_regions,cmap=plt_para['bnw_regions'],vmin=2,vmax=3,shading='auto')
    pimplausible = ax.pcolormesh(x_arr,y_arr,plot_implausible,cmap=plt_para['bnw_regions'],vmin=1,vmax=3,alpha=0.4,shading='auto')

    # if not plt_para['multi']:
    #
    #     #ax[i/2,i%2].pcolormesh(simulation[para_order[ax_list[1]]],simulation[para_order[ax_list[0]]],plot_no_peak,cmap=bnw,vmin=-1,vmax=2)
    #     #ax[i/2,i%2].pcolormesh(simulation[para_order[ax_list[1]]],simulation[para_order[ax_list[0]]],plot_inconsistent,cmap=bnw,vmin=-1,vmax=2)
    #
    #     mask_DM = ~np.isnan(results['DM_trans'])
    #     mask_no_peak = ~np.isnan(results['no_peak_trans'])
    #     mask_inc = ~np.isnan(results['inc_trans'])
    #
    #     ax.plot(results['DM_trans'][mask_DM],y_arr[mask_DM],'r-',linewidth=2,label=r'$\displaystyle \bar{\nu}_{DM}$')
    #     ax.plot(results['no_peak_trans'][mask_no_peak],y_arr[mask_no_peak],'k--',linewidth=2,label=r'$\displaystyle \bar{\nu}_{no\,peak}$')
    #     ax.plot(results['inc_trans'][mask_inc],y_arr[mask_inc],'k',linewidth=2,label=r'$\displaystyle \bar{\nu}_{inc}$')
    #     ax.plot(results['nu_implausible'][mask_inc],y_arr[mask_inc],'k:',linewidth=2,label=r'$\displaystyle \bar{\nu}_{imp}$')
    #
    #     # ax.set_xlabel(ax_labels[1],fontsize=12)
    #     # ax.set_ylabel(ax_labels[0],fontsize=12)

    ax.set_xlim([x_arr[0],x_arr[-1]])
    ax.set_ylim([y_arr[0],y_arr[-1]])
    return pchi,pgamma

def plot_colorbar(pchi,pgamma,plt_para,x=[0.85,0.88],y=[0.1,0.95]):

    if not plt_para['cb_plotted']:
        plt_para['cb_plotted'] ^= True;

        xL = x[1]-x[0]
        yL = y[1]-y[0]

        axcb1 = plt.axes([x[0],y[0]+yL*0.35,xL,yL*0.65])
        axcb2 = plt.axes([x[0],y[0],xL,yL*0.25])

        axcb1.tick_params(axis='both', which='major', labelsize=12)
        axcb2.tick_params(axis='both', which='major', labelsize=12)

        plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,4))
        plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,10),ticks=np.linspace(0,1,2))

        #axcb2.set_yticks(np.linspace(1,0,3))
        #axcb2.set_yticklabels(np.linspace(1,0,3))
        axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=12)
        axcb2.set_ylabel(r'$\displaystyle \gamma^2$',fontsize=12)

        plt.subplots_adjust(left=0.12, bottom=0.1, right=0.8, top=0.95, wspace=0.35, hspace=0.3)

def remove_frame(ax,positions=None):

    if positions is None:
        positions = ['left','right','top','bottom']

    for p in positions:
        ax.spines[p].set_visible(False)

def set_title(ax,order=1,title='',offset=(-0.1,0.1),pad=0):

    # pos_box = ax.get_position()
    # print(pos_box)
    # pos = [pos_box.x0,pos_box.y1]
    # # pos = fig.transFigure.transform(plt.get(ax,'position'))
    # x = pos[0]+offset[0]
    # y = pos[1]+offset[1]
    #
    # print(pos)
    ax.set_title(r'%s) %s'%(chr(96+order),title),position=offset,pad=pad)#,loc='left')#,fontsize=12)
    # ax.text(x=x,y=y,s=r'%s) %s'%(chr(96+order),title),ha='center',va='center',transform=ax.transAxes)#,loc='left')#,fontsize=12)
