import numpy as np

import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

from general.network import Network  # , distribution

def set_plot_params():
    # plot parameters
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    # plt.rcParams['font.family'] = ['Tahoma','Verdana']
    plt.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10


def plot_q(ax,res,plt_para,bound='imp',idxs=None,approx=False,order=1):

    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape

    if not idxs:
        idxs = range(steps1)

    ymax = 0
    ## plotting q
    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        trans_idx = get_trans_idx(res,bound,a,0,0)
        ymax = max(res['q'][0,a,:].max(),ymax)
        if trans_idx:
            ax.plot(res[x_key],res['q'][0,a,:],color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx],res['q'][0,a,:trans_idx],color=(col,0,0),label=r'$\alpha_0 = %g$'%res['alpha_0'][a])
        else:
            ax.plot(res[x_key],res['q'][0,a,],color=(col,0,0),label=r'$\alpha_0 = %g$'%res['alpha_0'][a])

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.1*steps)],res[x_key][:int(0.1*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (steps1 > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    ymax = 200
    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,ymax*1.1],
        # ylim=[0,20],
        xlabel=r'$\bar{\nu}\,$[Hz]',
        ylabel=r'$q\,$[Hz$^2$]'
    )

    if order:
        set_title(ax,order=order,title='second moment')


def plot_q_zoom(ax,res,plt_para,bound='imp',idxs=None,approx=False,order=1):

    x_key = plt_para['x']['key']
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
        ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        ax.plot(res[x_key][:int(0.012*steps)],res['rateWnt'][:int(0.012*steps)]*res['rate_max'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (len(res[x_key]) > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,5],
        xlabel=r'$\bar{\nu}\,$[Hz]',
        ylabel=r'$q\,$[Hz$^2$]'
    )

    if order:
        set_title(ax,order=order,title=r'the low $\bar{\nu}$ limit')


def plot_currents(ax,res,plt_para,bound='imp',idxs=None,approx=False,plot_options=['var','I'],order=1):
    
    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape
    # print(steps1,steps)
            
    ## plot constant lines
    # ax.axhline(0,c='k',ls='--',lw=0.5)

    if not idxs:
        idxs = range(steps1)

    idx_2 = np.argmin(np.abs(plt_para['x']['lim']/2-res[x_key]))
    idx_3 = np.argmin(np.abs(plt_para['x']['lim']/3-res[x_key]))
    if 'var' in plot_options:
    
        ## plot dark-matter transitions for one of the solutions
        trans_idx = get_trans_idx(res,'DM',2,0,0)
        if trans_idx:
            val = res[x_key][trans_idx]
            ax.axvline(val,color='k',ls=':')
            ax.text(val+0.2,0.2,r'$\bar{\nu}_{DM}$',fontsize=10)
        
        ## plotting temporal fluctuations
        trans_idx = get_trans_idx(res,bound,0,0,0)
        ax.plot(res[x_key][:trans_idx],res['sigma_V'][0,0,:trans_idx],color='k',ls='--',label=r'$\sigma_V$')
        ax.plot(res[x_key][trans_idx:],res['sigma_V'][0,0,trans_idx:],color='k',ls=':')
        
        ## plotting quenched fluctuations
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            trans_idx = get_trans_idx(res,bound,i,0,0)
            print(trans_idx)
            
            ax.plot(res[x_key][:trans_idx],res['alpha'][0,i,:trans_idx],color=(col,0,0),ls='-')#,label=r'$\alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
            ax.plot(res[x_key][trans_idx:],res['alpha'][0,i,trans_idx:],color=(col,0,0),ls=':')

        ax.text(res[x_key][idx_3],res['sigma_V'][0,0,idx_3]+0.05,r'$\sigma_{V_k}$',fontsize=10)
        ax.text(res[x_key][idx_2],res['alpha'][0,0,idx_2]-0.05,r'$\alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)

        plt.setp(ax, ylim=[-0.02,0.15])
        
    if 'I' in plot_options:
        
        ## plotting currents
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            trans_idx = get_trans_idx(res,bound,i,0,0)
            
            ax.plot(res[x_key],-res['I_balance'][0,i,:],':',color=[0.7,0,0],lw=0.8)
            ax.plot(res[x_key][:trans_idx],-res['I_balance'][0,i,:trans_idx],'-',color=(col,0,0),label=r'$I_{balance}$',lw=0.8)

        ax.text(res[x_key][idx_3],-res['I_balance'][0,0,idx_3]+0.05,r'$\bar{I}_0-\Psi_0$',fontsize=10)

        plt.setp(ax,
            ylim=[-0.3,0.02]
        )

    if 'var' in plot_options and 'I' in plot_options:
        plt.setp(ax,
            ylim=[-0.3,0.2]
        )

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        xlabel=r'$\bar{\nu}\,$[Hz]',
        ylabel=r'$I$'
    )

    if order:
        set_title(ax,order=order,title='variances')


def plot_gamma(ax,res,plt_para,bound='imp',idxs=None,approx=False,order=1):

    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape

    ax.plot(res[x_key],np.ones(steps),'k:',linewidth=0.5)

    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        # if a!=1:
        for p in range(Npop):   
            col = i/float(len(idxs)-1) if len(idxs)>1 else (Npop-p-1)/(Npop-1)
            trans_idx = get_trans_idx(res,bound,a,0,0)
            ax.plot(res[x_key],res['gamma'][p,a,:]**2,color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx],res['gamma'][p,a,:trans_idx]**2,color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['gamma_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    if order:
        set_title(ax,order=order,title=r'$\quad$ DM exponent $\gamma$')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,7],
        xlabel=get_displayString(plt_para['x']['key']),
        ylabel=r'$\gamma^2$'
    )


def plot_chi(ax,res,plt_para,bound='imp',idxs=None,approx=False,order=1):

    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape

    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        # p=0
        for p in range(Npop):
            col = i/float(len(idxs)-1) if len(idxs)>1 else (Npop-p-1)/(Npop-1)

            trans_idx = get_trans_idx(res,bound,a,0,p)
            ax.plot(res[x_key],res['chi'][p,a,:],color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx],res['chi'][p,a,:trans_idx],color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['chi_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,10],
        xlabel=get_displayString(plt_para['x']['key']),
        ylabel=r'$\chi$'
    )

    if order:
        set_title(ax,order=order,title=r'skewness coeff. $\chi$')


def plot_regions(ax,res,plt_para,order=1):

    # Npop,steps1,steps = res['q'].shape

    # mask_DM = get_trans_idx(res,'DM',0,0,0)
    # print(mask_DM)
    mask_DM = ~res['DM_trans'][0,:,0].mask
    # print(mask_DM)
    mask_no_peak = ~res['np_trans'][0,:,0].mask
    mask_inc = ~res['inc_trans'][:,0].mask if np.any(res['inc_trans'][:,0].mask) else np.ones_like(res['inc_trans'][:,0],'bool')
    mask_imp = ~res['imp_trans'][:,0].mask if np.any(res['imp_trans'][:,0].mask) else np.ones_like(res['imp_trans'][:,0],'bool')

    x_arr = res[plt_para['x']['key']]

    # print(x_arr[res['DM_trans'][0,mask_DM,0]])
    # print(res['alpha_0'][mask_DM])

    ax.plot(x_arr[res['DM_trans'][0,mask_DM,0]],res['alpha_0'][mask_DM],'r-',label=r'$\bar{\nu}_{DM}$')
    ax.plot(x_arr[res['np_trans'][0,mask_no_peak,0]],res['alpha_0'][mask_no_peak],'k--',label=r'$\bar{\nu}_{no\,peak}$')
    ax.plot(x_arr[res['inc_trans'][mask_inc,0]],res['alpha_0'][mask_inc],'k',label=r'$\bar{\nu}_{inc}$')
    ax.plot(x_arr[res['imp_trans'][mask_imp,0]],res['alpha_0'][mask_imp],'k:',label=r'$\bar{\nu}_{imp}$')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,0.2],
        xlabel=r'$\bar{\nu}\,$[Hz]',
        ylabel=r'$\alpha_0$'
    )

    # ax.set_yticks(np.linspace(0,0.2,6))

    ax.legend(prop={'size':10},ncol=2,bbox_to_anchor=(0.05,1.4),loc='upper left',handlelength=1)

    if order:
        set_title(ax,order=order,title='boundaries')


def plot_fins(ax,x_arr,y_arr,gamma,chi,regions,implausible=False):

    levs = 20
    plt_para = {}
    plt_para['cb_plotted'] = False
    plt_para['bnw'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=levs-1)
    plt_para['heat'] = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=levs-1)
    plt_para['bnw_regions'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)

    mask_inconsistent = (regions == 3)
    mask_no_peak = (regions == 2)
    if np.any(implausible==1):
        mask_implausible = (implausible == 1)
    else:
        mask_implausible = np.zeros_like(regions,'bool')

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
        # ax[i/2,i%2].pcolormesh(simulation[para_order[ax_list[1]]],simulation[para_order[ax_list[0]]],plot_inconsistent,cmap=bnw,vmin=-1,vmax=2)
        # ax.pcolormesh(x_arr,y_arr,plot_inconsistent,cmap=plt_para['bnw'],vmin=-1,vmax=2)
    #
    #     mask_DM = ~np.isnan(results['DM_trans'])
    #     mask_no_peak = ~np.isnan(results['no_peak_trans'])
    #     mask_inc = ~np.isnan(results['inc_trans'])
    #
    #     ax.plot(results['DM_trans'][mask_DM],y_arr[mask_DM],'r-',linewidth=2,label=r'$\bar{\nu}_{DM}$')
    #     ax.plot(results['no_peak_trans'][mask_no_peak],y_arr[mask_no_peak],'k--',linewidth=2,label=r'$\bar{\nu}_{no\,peak}$')
    #     ax.plot(results['inc_trans'][mask_inc],y_arr[mask_inc],'k',linewidth=2,label=r'$\bar{\nu}_{inc}$')
    #     ax.plot(results['nu_implausible'][mask_inc],y_arr[mask_inc],'k:',linewidth=2,label=r'$\bar{\nu}_{imp}$')
    #
    #     # ax.set_xlabel(ax_labels[1],fontsize=12)
    #     # ax.set_ylabel(ax_labels[0],fontsize=12)
    # x_arr[0] = 0
    # y_arr[0] = 0
    ax.set_xlim([x_arr[0],x_arr[-1]])
    ax.set_ylim([y_arr[0],y_arr[-1]])
    return pchi,pgamma


def plot_3D(ax,results,para_order,sim_vals):
    
    ## initialize parameters and arrays
    px,py,pz = para_order
    
    nP,steps,_ = results[0]['gamma'].shape
    sim_steps = len(results)
    
    data_3D_inc = np.full([sim_steps,steps,3],np.nan)
    data_3D_DM = np.full([sim_steps,steps,3],np.nan)
    
    cross_point = np.zeros((sim_steps,3))

    p=0
    
    ## iterate through each of the results
    for i,res in enumerate(results):
        
        ## define dark matter transition lines
        DM_mask = ~res['DM_trans'][p,:,0].mask
        # data_mask = ~res['DM_trans'][p,:,0].mask
        N_DM_entries = DM_mask.sum()
        DM_idxes = res['DM_trans'][p,DM_mask,0]
        other_idxes = np.where(DM_mask)[0]
        
        data_3D_DM[i,:N_DM_entries,0] = res[px][DM_idxes]
        data_3D_DM[i,:N_DM_entries,1] = sim_vals[i]
        data_3D_DM[i,:N_DM_entries,2] = res[pz][other_idxes]
        # data_3D_DM[i,:N_DM_entries,1] = res[py][other_idxes] if py in order else sim_vals[i]
        # data_3D_DM[i,:N_DM_entries,2] = res[pz][other_idxes] if pz in order else sim_vals[i]

        ## define inconsistency transition lines
        inc_mask = ~res['inc_trans'][:,0].mask if np.any(res['inc_trans'].mask) else np.ones((steps),'bool')
        N_inc_entries = inc_mask.sum()
        inc_idxes = res['inc_trans'][inc_mask,0]
        other_idxes = np.where(inc_mask)[0]
        
        data_3D_inc[i,:N_inc_entries,0] = res[px][inc_idxes]
        data_3D_inc[i,:N_inc_entries,1] = sim_vals[i]
        data_3D_inc[i,:N_inc_entries,2] = res[pz][other_idxes]
        # data_3D_inc[i,:N_inc_entries,1] = res[py][other_idxes] if py in order else sim_vals[i]
        # data_3D_inc[i,:N_inc_entries,2] = res[pz][other_idxes] if pz in order else sim_vals[i]
        
        
        ## remove data above DM transition to maintain clean plots
        cross_point[i,:] = data_3D_DM[i,:N_DM_entries,:][-1,:]
        # data_3D_inc[i,data_3D_inc[i,:,2]>cross_point[i,2],:] = np.nan
        
        ## plot transition lines of DM and inc (only plot every 2nd data)
        if i%2==0 and sim_steps>5:
            ax.plot(data_3D_DM[i,:,0],data_3D_DM[i,:,1],data_3D_DM[i,:,2],'k-',lw=1)
            ax.plot(data_3D_inc[i,:,0],data_3D_inc[i,:,1],data_3D_inc[i,:,2],'k-',lw=1)
        else:
            ax.plot(data_3D_DM[i,:,0],data_3D_DM[i,:,1],data_3D_DM[i,:,2],'k:',lw=0.5)
            ax.plot(data_3D_inc[i,:,0],data_3D_inc[i,:,1],data_3D_inc[i,:,2],'k:',lw=0.5)
    
    ## plot connecting lines between different slices
    ax.plot(data_3D_DM[:,0,0],data_3D_DM[:,0,1],data_3D_DM[:,0,2],'k-',lw=1)
    ax.plot(data_3D_inc[:,0,0],data_3D_inc[:,0,1],data_3D_inc[:,0,2],'k-',lw=1)
    ax.plot(cross_point[:,0],cross_point[:,1],cross_point[:,2],'k-')
    
    ## plot surface area of DM and inc transitions
    for i in range(sim_steps-1):
        
        ## DM
        X = np.concatenate((data_3D_DM[i,:,0],np.flipud(data_3D_DM[i+1,:,0])))
        mask = ~np.isnan(X)
        if np.any(mask):
            Y = np.concatenate((data_3D_DM[i,:,1],np.flipud(data_3D_DM[i+1,:,1])))
            Z = np.concatenate((data_3D_DM[i,:,2],np.flipud(data_3D_DM[i+1,:,2])))
            X = np.pad(X[mask],1,mode='wrap')
            Y = np.pad(Y[mask],1,mode='wrap')
            Z = np.pad(Z[mask],1,mode='wrap')
            ax.plot_trisurf(X,Y,Z,color='red',alpha=0.3)
        
        ## inc
        X = np.concatenate((data_3D_inc[i,:,0],np.flipud(data_3D_inc[i+1,:,0])))
        mask = ~np.isnan(X)
        if np.any(mask):
            Y = np.concatenate((data_3D_inc[i,:,1],np.flipud(data_3D_inc[i+1,:,1])))
            Z = np.concatenate((data_3D_inc[i,:,2],np.flipud(data_3D_inc[i+1,:,2])))
            X = np.pad(X[mask],1,mode='wrap')
            Y = np.pad(Y[mask],1,mode='wrap')
            Z = np.pad(Z[mask],1,mode='wrap')
            # Z = np.pad(Z[mask],1,mode='constant',constant_values=0)
            ax.plot_trisurf(X,Y,Z,color='grey',alpha=0.3)
    
    max_x_DM, max_y_DM, max_z_DM = np.nanmax(data_3D_DM,axis=(0,1))
    max_x_inc, max_y_inc, max_z_inc = np.nanmax(data_3D_inc,axis=(0,1))
    
    max_x = max(max_x_DM,max_x_inc)
    max_y = max(max_y_DM,max_y_inc)
    max_z = max(max_z_DM,max_z_inc)
    plt.setp(ax,
         xticks=np.linspace(0,20,5),
         yticks=np.linspace(0,0.1,3),yticklabels=np.linspace(0,100,3),
         zticks=np.linspace(0,0.1,3),
         xlim=[0,max_x],
         ylim=[0,max_y],
         zlim=[0,max_z],
    )


def plot_colorbar(pchi,pgamma,x=[0.95,0.97],y=[0.1,0.95]):

    # if not plt_para['cb_plotted']:
        # plt_para['cb_plotted'] ^= True

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
    axcb1.set_ylabel(r'$\chi$',fontsize=12)
    axcb2.set_ylabel(r'$\gamma^2$',fontsize=12)

    #plt.subplots_adjust(left=0.12, bottom=0.1, right=0.8, top=0.95, wspace=0.35, hspace=0.3)


def plot_approx(ax_ex,ax_I,ax_alpha):

    res_approx = get_approximation(0.00)
    res_approx_hetero = get_approximation(0.04)
    # print(res['q_approx'][0,0,:])
    # print(res_approx['q'])
    #
    # print(res['q_approx'][0,2,:])
    # print(res_approx_hetero['q'])

    net = network(alpha_0=0.00)
    ylims = [-0.25,0]
    nu = 5
    p=0
    q = np.linspace(0,4*nu**2,1001)
    y_is = -np.sqrt(net.I_squared_nu(nu,nu**2,p)) ## intersection with function
    y_is_transformed = (y_is-ylims[0]) / (ylims[1] - ylims[0])
    min_idx = np.nanargmin(abs(np.sqrt(net.I_squared_q(nu,q,p))-np.sqrt(net.I_squared_nu(nu,q,p))))

    ax_ex.axvline(nu**2,ymax=y_is_transformed,color=[0.6,0.6,0.6],lw=2,ls='--')
    ax_ex.axhline(-np.sqrt(net.I_squared_nu(nu,nu**2,p)),color=[0.6,0.6,0.6],lw=2,ls='--')#,label=r"approx.: $q=\bar{\nu}^2$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_nu(nu,q,p)),'k-',label=r"solution for $\bar{\nu}$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_q(nu,q,p)),'k--',label=r"solution for $q$")
    ax_ex.annotate(r'$q = \bar{\nu}^2$',xy=[nu**2,-np.sqrt(net.I_squared_nu(nu,nu**2,p))],xytext=[nu**2-22,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)
    ax_ex.annotate(r'$(\bar{\nu}^{\star},q^{\star})$',xy=[q[min_idx],-np.sqrt(net.I_squared_nu(nu,q[min_idx],p))],xytext=[nu**2+15,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)

    plt.setp(ax_ex, ylim=ylims)
    ax_ex.legend(prop={'size':10},bbox_to_anchor=(1.2,0.0),loc='lower right',frameon=False)
    loc = ticker.MultipleLocator(base=0.1) # this locator puts ticks at regular intervals
    ax_ex.yaxis.set_major_locator(loc)
    plot_I(ax_I,res_approx,'k')
    plot_I(ax_I,res_approx_hetero,'r')

    plt.setp(ax_I,xlabel=r'$\bar{\nu}$ [Hz]',ylabel=r'$\frac{|\Delta I|}{\Psi_0 - I_0}$')

    ax_alpha.plot(res_approx['nu'],res_approx['alpha'][:,0],'k')
    ax_alpha.plot(res_approx['nu'],res_approx['alpha'][:,1],'k--')
    ax_alpha.plot(res_approx['nu'],res_approx_hetero['alpha'][:,0],'r')
    ax_alpha.plot(res_approx['nu'],res_approx_hetero['alpha'][:,1],'r--')


def get_approximation(alpha=0.0):
    net = network(eps=0,tau_G=0.01,rateWnt=5,alpha_0=alpha)

    nu_lim = 12
    steps = 101
    # I_exact = np.zeros(steps)
    res = {
        'nu': np.linspace(1/steps,nu_lim,steps),
        'I': np.zeros((steps,2)),
        'q': np.zeros((steps,2)),
        'alpha': np.zeros((steps,2)),
    }

    for i,nu in enumerate(res['nu']):


        q = np.linspace(0,max(10,4*nu**2),1001)
        dI = abs(np.sqrt(net.I_squared_q(nu,q,0))-np.sqrt(net.I_squared_nu(nu,q,0)))
        min_idx = np.nanargmin(dI[1:])
        min_val = np.nanmin(dI[1:])
        #print(dI)#
        #print(min_idx)

        # print(i, nu, min_idx, min_val)
        #print(nu,nu**2)
        #print(net.get_q(nu,nu**2,0,I=np.sqrt(net.I_squared_nu(nu,q[min_idx],0))))
        res['I'][i,0] = np.sqrt(net.I_squared_nu(nu,q[min_idx],0))
        res['I'][i,1] = np.sqrt(net.I_squared_nu(nu,nu**2,0))
        res['q'][i,0] = q[min_idx]
        res['q'][i,1:] = net.get_q(nu,nu**2,0,I=res['I'][i,0]),
        res['alpha'][i,0] = net.alpha(res['q'][i,0],0)
        res['alpha'][i,1] = net.alpha(res['q'][i,1],0)


    # fig,axes = plt.subplots(2,2)
    # axes[0,0].plot(res['q'][:,0],'k-')
    # axes[0,0].plot(res['q'][:,1],'k--')
    #
    # axes[1,0].plot(res['alpha'][:,0],'k')
    # axes[1,0].plot(res['alpha'][:,1],'k--')
    # plt.setp(axes[1,0],ylim=[0,0.15])
    #
    # axes[0,1].plot(res['I'][:,0],'k')
    # axes[0,1].plot(res['I'][:,1],'k--')
    # plt.show(block=False)

    return res

def plot_I(ax,res,ls):
    # nu_lim = 12
    # steps = res['I'].shape[0]
    ax.plot(res['nu'],np.abs(res['I'][:,0]-res['I'][:,1])/res['I'][:,0],ls,label=r"$\frac{\Delta I}{I}$")

    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_exact,'k')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,0],'r')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),q_approx[:,1],'r--')
    # ax_inset.plot(np.linspace(0,nu_lim,steps),I_approx,'r')


# def plot_distr(alpha_0=[0,0.04],rateWnt=[0.5,2,10],tau_M=0.010,save=0,file_format='png'):
def plot_distr(ax,idx_rate,res,plt_para,plt_style="row",annotate=['highlight_var'],show_xticks=True,order=1):
    '''
        contains plotting routing for a single solution to the network firing rate distribution
        shows...
         ... input current distribution
         ... f-I curve / firing rate response
         ... firing rate distribution
        
        allows for plotting as "row" or "sketch" plt_style
        annotate:
            0: no annotations
            1: highlight variances without labels
            2: highlight variances and add labels
    '''

    steps = 10000
    p = 0
    rateWnt = res['rateWnt'][idx_rate]

    ### define axes

    ## get location from subplots
    loc = ax.get_position().get_points()
    ax.remove()

    ## define subaxes
    if plt_style=='row':
        y_offset = loc[0,1]
        y_height = loc[1,1]-loc[0,1]

        ax_fI = plt.axes(position=[0.1,y_offset,0.65,y_height])
        ax_I = ax_fI.twinx()
        ax_distr = plt.axes(position=[0.8,y_offset,0.15,y_height])

        I_range = np.linspace(-0.35,0.05,steps+1)

        ## set plot parameters
        for axx in [ax_fI,ax_I,ax_distr]:
            axx.spines[['top','right']].set_color('none')
            axx.yaxis.set_ticks_position('left')
            axx.xaxis.set_ticks_position('bottom')

    else:
        x_offset = [0.02,0.3,0.98]
        y_offset = [0.1,0.5,0.925]

        x_width_right = x_offset[2] - x_offset[1]
        x_width_left = x_offset[1] - x_offset[0] #- 0.05
        y_height_top = y_offset[2] - y_offset[1]
        y_height_bottom = y_offset[1] - y_offset[0]

        ax_I = plt.axes(position=[x_offset[1],y_offset[0],x_width_right,y_height_bottom])
        ax_fI = plt.axes(position=[x_offset[1],y_offset[1],x_width_right,y_height_top],sharex=ax_I)
        ax_distr = plt.axes(position=[x_offset[0],y_offset[1],x_width_left,y_height_top])

        ax_distr.spines[['top','left']].set_color('none')
        # ax_distr.yaxis.set_ticks_position('right')
        plt.setp(ax_distr,yticks=[])

        for axx in [ax_fI,ax_I]:
            axx.spines[['top','left','right']].set_color('none')
            plt.setp(axx,xticks=[])

        I_range = np.linspace(-0.45,0.3,steps+1)

        marker_lineprops = {'color':'k','linewidth':0.7,'linestyle':'--'}
        I_marked_position = 0.16
        I_marked_interval = 0.015

        I_marked = np.zeros((2,2))
        for j,jj in enumerate([-1,1]):
            for i,ii in enumerate([-1,1]):
                I_marked[j,i] = jj*I_marked_position-jj*ii*I_marked_interval

    ax_I.set_yticks([])

    ## set some general parameters
    factor = 0.6
    I_max = 1.

    if not show_xticks:
        plt.setp(ax_fI,xticklabels=[])

    nu_range_max = 0

    ## plot f-I curve / firing rate response
    idxes = (p,0,idx_rate)
    nu_max = res['rate_max'][idxes]
    fI = nu_max * np.exp(-I_range**2/(2*res['sigma_V'][idxes]**2))
    ax_fI.plot(I_range,fI,'k')

    ax_fI.plot([0,0],[0,nu_max],'k--',linewidth=0.5)

    if 'highlight_var' in annotate:
        ax_fI.plot([-res['sigma_V'][idxes],0],[factor*nu_max,factor*nu_max],'k',linewidth=2)

        if 'highlight_var_text' in annotate:
            ax_fI.text(0.005,1,r'$\Psi_0$',fontsize=12)
            ax_fI.text(-res['sigma_V'][idxes]/2.-0.01,factor*nu_max-2,r'$\sigma_V$',fontsize=10)

    idxes_alpha = [0,2] if plt_style=='row' else [0]

    for idx_alpha, alpha_0 in enumerate(res['alpha_0'][idxes_alpha]):

        idxes = (p,idx_alpha,idx_rate)

        if idx_alpha==0:
            lineprops = {'linewidth':1.5,'linestyle':'-'}
            I_col = 'grey'
        elif idx_alpha==1:
            lineprops = {'linewidth':1.5,'linestyle':':'}
            I_col = 'r'

        ## plot input current distribution
        I_balance = -res['I_balance'][idxes]
        p_I = np.exp(-(I_range-I_balance)**2/(2*res['alpha'][idxes]**2))
        ax_I.plot(I_range,p_I,color=I_col,**lineprops)

        ax_I.plot([I_balance,I_balance],[0,I_max],I_col,linewidth=1.,linestyle=':')
        if 'highlight_var' in annotate:
            ax_I.plot([I_balance-res['alpha'][idxes],I_balance],[factor*I_max,factor*I_max],I_col,**lineprops,label=r'$\alpha_I$')

            if ('highlight_var_text' in annotate) and (idx_alpha==0):
                ax_I.text(-res['I_balance'][idxes]-res['alpha'][idxes]/2.-0.01,factor*I_max - 0.15,r'$\alpha_I$',fontsize=10)
                ax_I.text(-res['I_balance'][idxes]+0.005,0.05,r'$\bar{I}_0$',fontsize=10)

        ## plot firing rate distribution
        range_rate = np.linspace(0,nu_max,steps+1)
        p_nu = distribution(range_rate,res['gamma'][idxes],res['delta'][idxes],res['rate_max'][idxes])
        rho_at_mean = p_nu[np.argmin(abs(range_rate-rateWnt))]

        if plt_style=='row':
            ax_distr.plot([0,rho_at_mean],[rateWnt,rateWnt],'k--')

        p_nu_max = np.nanmax(p_nu)
        nu_range_max = max(nu_range_max,range_rate[np.where(p_nu>10**(-3)*p_nu_max)[0][-1]])*1.1

        ax_distr.plot(p_nu,range_rate,linewidth=1,linestyle='-',color='k' if idx_alpha==0 else 'r')

        plt.setp(ax_distr,xlim=[0,p_nu_max*1.05])

    if plt_style=='sketch':

        idxes_rate = np.zeros(2,'int')
        for j,I_marked_row in enumerate(I_marked):

            for i,I_mark in enumerate(I_marked_row):
                ax_I.axvline(I_mark,**marker_lineprops)

                nu_marked = fI[np.argmin(abs(I_range-I_mark))]
                ax_fI.plot([I_mark,I_mark],[0,nu_marked],**marker_lineprops)

                if j==0:
                    print(nu_marked)
                    ax_I.text(I_mark-0.015+(i*0.02),-0.15,f'$I_{i+1}$',fontsize=8)
                    ax_fI.text(I_range[0]+0.01,nu_marked-2.5+(i*4),'$\\nu(I_{%i})$'%(i+1),fontsize=8)
                else:
                    ax_fI.plot([I_range[0],I_mark],[nu_marked,nu_marked],**marker_lineprops)

                    idxes_rate[i] = np.argmin(abs(range_rate-nu_marked))

                    p_nu_marked = p_nu[idxes_rate[i]]
                    ax_distr.plot([p_nu_marked,0],[nu_marked,nu_marked],**marker_lineprops)

            idxes_I = [np.argmin(abs(I_range-I_marked_row[0])),np.argmin(abs(I_range-I_marked_row[1]))]
            idxes_I.sort()# = [min(idxes_I),max(idxes_I)]

            ax_I.fill_between(I_range[idxes_I[0]:idxes_I[1]],0,p_I[idxes_I[0]:idxes_I[1]],color='grey',alpha=0.6)

            ax_fI.fill_between(I_range[idxes_I[0]:idxes_I[1]],0,fI[idxes_I[0]:idxes_I[1]],color='grey',alpha=0.6)

        idxes_rate.sort()
        ax_distr.fill_betweenx(range_rate[idxes_rate[0]:idxes_rate[1]],0,p_nu[idxes_rate[0]:idxes_rate[1]],color='grey',alpha=0.6)

        ax_fI.plot([I_range[0],0],[nu_max,nu_max],linewidth=1.,linestyle=':',color='k')
        ax_distr.plot([p_nu_max/2.,0],[nu_max,nu_max],linewidth=1.,linestyle=':',color='k')

        ax_distr.text(p_nu_max/2.,nu_max+1,'$\\bar{\\nu}^{max}$',fontsize=10)
        ax_distr.annotate('pole',xy=[0,nu_max],xytext=[p_nu_max/4.,nu_max-5],arrowprops=dict(arrowstyle="->"),fontsize=10)

        ax_I.set_xlabel('$I$', loc="right")
        ax_fI.set_ylabel('$\\nu$', loc="top", rotation=0)
        ax_distr.set_xlabel('$\\rho(\\nu)$', loc="left")

        ax_fI.yaxis.set_label_coords(0.0, 1.02)

    plt.setp(ax_I,ylim=[0,1.1])
    plt.setp(ax_distr,xticks=[],yticklabels=[])

    if plt_style=='row':
        for axx in [ax_fI,ax_distr]:
            plt.setp(axx,yticks=np.linspace(0,(nu_max//5)*5,5),ylim=[0,nu_range_max])

        ax_distr.text(0.3*p_nu_max,nu_range_max*3/4.,r'$\chi \approx %4.2f \rightarrow %4.2f $'%(res['chi'][p,0,idx_rate],res['chi'][p,1,idx_rate]),bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=10)

        plt.setp(ax_fI,xlim=I_range[[0,-1]],xticks=np.linspace(-0.3,0.0,4))
    else:
        for axx in [ax_fI,ax_distr]:
            plt.setp(axx,ylim=[0,nu_max*1.05],yticks=[])

        plt.setp(ax_fI,xlim=I_range[[0,-1]])

        xlim = ax_distr.get_xlim()
        plt.setp(ax_distr,xlim=[xlim[1],0])

    if 'title' in annotate:
        ax_fI.text(I_range[int(steps*0.05)],nu_range_max*0.85,r'$\bar{\nu}=%g\,$Hz'%rateWnt,fontsize=12)

    return 

    # title_x = -0.1
    # title_y = 1.05

    nu_border = [5,15,25]

    for i in range(len(rateWnt)):
        for j in range(len(alpha_0)):
            results = get_samples_from_theory(tau_M=tau_M,T=1000,dftheory=0,dftime=0,p_theory=2,plot=1,rate=rateWnt[i],alpha_0=alpha_0[j])
            # print results
            # ax1 = plt.subplot(gs[2*i,0])
            if (j == 0):
                border_y = (len(rateWnt)-i-1)*v_plots + (len(rateWnt)-i-1)*v_spaces + v_bottom
                # ax1 = plt.axes([0.1,border_y+(v_plots+v_spaces_small)/2.,0.33,(v_plots-v_spaces_small)/2.])
                # ax2 = plt.axes([0.1,border_y,0.33,(v_plots-v_spaces_small)/2.])
                ax2 = plt.axes([0.1,border_y,0.45,v_plots])
                ax1 = ax2.twinx()
                ax3 = plt.axes([0.6,border_y,0.15,v_plots])

                ax2.set_xticks(np.linspace(-0.4,0.0,5))
                ax2.set_yticks([])
                ##ax2.set_xlabel(r'$I$')
                # ax2.set_ylabel(r'$\rho(I)$')
                ax2.spines['right'].set_color('none')
                ax2.yaxis.set_ticks_position('left')
                ax2.spines['top'].set_color('none')
                ax2.xaxis.set_ticks_position('bottom')

                # ax1.set_xticks([])
                # ax1.set_ylabel(r'$\nu(I)\,$[Hz]')
                ax1.spines['right'].set_color('none')
                ax1.yaxis.set_ticks_position('left')
                ax1.spines['top'].set_color('none')
                ax1.xaxis.set_ticks_position('bottom')

                ax3.set_xticks([])
                ax3.set_yticks([])
                # ax3.set_xlabel(r'$\nu\,$[Hz]')
                # ax3.set_ylabel(r'$\rho(\nu)$')
                ax3.spines['right'].set_color('none')
                ax3.yaxis.set_ticks_position('left')
                ax3.spines['top'].set_color('none')
                ax3.xaxis.set_ticks_position('bottom')

                # factor = 1-1./math.pi
                factor = 0.6

                nu_max = max(results['f_I'])
                print(nu_max)
                ax1.plot([-results['sigma_V'],0],[factor*nu_max,factor*nu_max],'k',linewidth=2)#,label=r'$\sigma_V$')
                ax1.plot(results['I_range'],results['f_I'],'k')
                ax1.plot([0,0],[0,nu_max],'k--',linewidth=0.5)

                ax1.set_yticks(np.linspace(0,25,6))

                # ax1.legend(prop={'size':10},loc=3)

                I_max = max(results['I_distr'])
                I_tmp = np.copy(I_max)
                ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r',linewidth=2,label=r'$\alpha_I$')
                ax2.plot(results['I_range'],results['I_distr'],'r')
                ax2.plot([results['I'],results['I']],[0,I_tmp],'r--',linewidth=0.5)
                # ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r',linewidth=2,label=r'$\alpha_I$')
                # ax2.plot(results['I_range'],results['I_distr'],'k')
                # ax2.plot([results['I'],results['I']],[0,I_max],'k--',linewidth=0.5)
                # ax2.plot([0,0],[0,1.1*I_max],'k--',linewidth=0.5)

                # ax2.legend(prop={'size':10},loc=4)

                # print results['p_range']
                rho_at_mean = results['p_exact'][np.argmin(abs(results['p_range']-rateWnt[i]))]
                rho_max1 = max(results['p_exact'])

                # ax3.plot([rateWnt[i],rateWnt[i]],[0,rho_at_mean],'k--')
                # ax3.plot(results['p_range'],results['p_exact'],'k',label='exact solution')
                # ax3.plot(results['p_range'],results['p_approx'],'r--',linewidth=1,label='approx. solution')
                # ax3.set_xlim([0,nu_border[i]])
                ax3.plot([0,rho_at_mean],[rateWnt[i],rateWnt[i]],'k--')
                ax3.plot(results['p_exact'],results['p_range'],'k',label='exact solution')
                ax3.plot(results['p_approx'],results['p_range'],'r--',linewidth=1,label='approx. solution')

                chi1 = np.copy(results['chi'])
                # print chi1
                if (i==0):
                    ax1.set_title(r'a)',position=(title_x,title_y))#,loc='left')
                    ax3.set_title(r'b)',position=(title_x,title_y))# homog.: $\alpha_0 = 0$',loc='left')
                    ax1.text(0.01,1,r'$\Psi$',fontsize=12)
                    y_border = 12
                else:
                    y_border = 25

                if (i==1):
                    # if results['sigma_V'] > 0.05:
                    ax1.text(-results['sigma_V']/2.-0.01,factor*nu_max*0.7,r'$\sigma_V$',fontsize=10)
                    # if results['alpha'] > 0.05:
                    ax2.text(results['I']-results['alpha']/2.-0.01,factor*I_max*0.7,r'$\alpha_I$',fontsize=10)

                if (i==2):
                    ax1.text(results['I']+0.01,1,r'$\bar{I}_0$',fontsize=10)

                ax1.set_ylim([0,y_border])
                ax2.set_ylim([0,1.1*I_max*y_border/25.])
                ax3.set_ylim([0,y_border])

                ax1.text(-0.45,0.9*y_border,r'$\bar{\nu}=%g\,$Hz'%rateWnt[i],fontsize=12)

                # ax3.text(rateWnt[i]+0.2,rho_at_mean,r'$\bar{\nu}$',fontsize=12)
                # else:
                # ax3.text(rho_at_mean,rateWnt[i]+0.05*nu_border[i],r'$\bar{\nu}$',fontsize=12)
            if (j == 1):
                # ax1.plot(results['I_range'],results['f_I'],'k')
                I_max = max(results['I_distr'])
                # ax2.plot(results['I_range'],results['I_distr']*I_tmp/I_max,'k:',linewidth=0.5)
                ax2.plot(results['I_range'],results['I_distr']*I_tmp/I_max,'r:',linewidth=0.5)
                ax2.plot([results['I'],results['I']],[0,I_tmp],'r:',linewidth=0.5)
                # factor = (1-1./np.exp(1))*0.9
                # ax2.plot([results['I']-results['alpha'],results['I']],[factor*I_max,factor*I_max],'r--',linewidth=2)
                ax4 = plt.axes([0.8,border_y,0.15,v_plots])

                ax4.set_xticks([])
                ax4.set_yticks([])
                # ax4.set_xlabel(r'$\nu\,$[Hz]')
                # ax4.set_ylabel(r'$\rho(\nu)$')
                ax4.spines['right'].set_color('none')
                ax4.yaxis.set_ticks_position('left')
                ax4.spines['top'].set_color('none')
                ax4.xaxis.set_ticks_position('bottom')

                rho_at_mean = results['p_approx'][np.argmin(abs(results['p_range']-rateWnt[i]))]
                rho_max2 = max(results['p_exact'])
                # print results['p_exact']
                if (rho_max2 == results['p_exact'][1]):
                    rho_max2 /= 3
                rho_max = max(rho_max1,rho_max2)
                # ax4.plot([rateWnt[i],rateWnt[i]],[0,rho_at_mean],'k--')
                # ax4.plot(results['p_range'],results['p_exact'],'k')
                # ax4.plot(results['p_range'],results['p_approx'],'r--')
                # ax4.set_xlim([0,nu_border[i]])
                ax4.plot([0,rho_at_mean],[rateWnt[i],rateWnt[i]],'k--')
                ax4.plot(results['p_exact'],results['p_range'],'k')
                ax4.plot(results['p_approx'],results['p_range'],'r--')

                if (i==0):
                    ax4.set_title(r'c)',position=(title_x,title_y))# inhom.: $\alpha_0 = 0.04$',loc='left')
                    ax4.legend(prop={'size':10},bbox_to_anchor=(-0.4,0.85),loc='lower left')

                ax4.set_ylim([0,y_border])

                ax3.text(0.6*rho_max,y_border*3/4.,r'$\chi \approx %4.2f$'%chi1,bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=10)
                ax4.text(0.6*rho_max,y_border*3/4.,r'$\chi \approx %4.2f$'%results['chi'],bbox={'facecolor':'white','alpha':0.9,'pad':5},fontsize=10)

                # else:
                # ax3.text(12.5,0.7*rho_max,r'$\chi \approx %4.2g$'%chi1,bbox={'facecolor':'white','alpha':0.9,'pad':5})
                # ax4.text(12.5,0.7*rho_max,r'$\chi \approx %4.2g$'%results['chi'],bbox={'facecolor':'white','alpha':0.9,'pad':5})
                ax3.set_xlim([0,rho_max])
                ax4.set_xlim([0,rho_max])

            ax2.set_xlim([-0.5,0.1])
            ax1.set_xlim([-0.5,0.1])
            # ax2.set_xlim([-0.5,0.2])
            ax1.set_ylim([0,y_border])

    plt.show(block=False)


def remove_frame(ax,positions=None):

    if positions is None:
        positions = ['left','right','top','bottom']

    for p in positions:
        ax.spines[p].set_visible(False)

def set_title(ax,order=1,title='',offset=(-0.1,1.1),pad=0,fontsize=10):

    # pos_box = ax.get_position()
    # print(pos_box)
    # pos = [pos_box.x0,pos_box.y1]
    # # pos = fig.transFigure.transform(plt.get(ax,'position'))
    # x = pos[0]+offset[0]
    # y = pos[1]+offset[1]
    #
    # print(pos)
    # ax.set_title(r'%s) %s'%(chr(96+order),title),position=offset,pad=pad,fontsize=fontsize)#,loc='left')#,fontsize=12)
    ax.set_title(r'%s) %s'%(chr(96+order),title),x=offset[0],y=offset[1],pad=pad,fontsize=fontsize)#,loc='left')#,fontsize=12)
    # ax.text(x=x,y=y,s=r'%s) %s'%(chr(96+order),title),ha='center',va='center',transform=ax.transAxes)#,loc='left')#,fontsize=12)


def get_displayString(key):

    if key == 'rateWnt':
        return "$\\bar{\\nu}\,$[Hz]"
    elif key == 'alpha_0':
        return r'$\alpha_0$'
    elif key == 'tau_I':
        return r'$\tau_I\,$[ms]'
    elif key == 'tau_G':
        return r'$\tau_G\,$[ms]'
    elif key == 'eps':
        return r'$\varepsilon$'
    elif key == 'eta':
        return r'$\eta$'
    elif key == 'tau_n':
        return r'$n$'
    elif key == 'Psi_0':
        return r'$\Psi_0$'
    else:
        return r''


def get_trans_idx(res,bound,y,n=0,p=0):

    if not ((bound+'_trans') in res.keys()):
        return None
    if bound in ['imp','inc']:
        mask = res[bound+'_trans'].mask
        mask_val = mask if np.prod(mask.shape)==1 else mask[y,n]
        
        # mask = np.all(res[bound+'_trans'].mask)
        return res[bound+'_trans'][y,n] if ~mask_val else None
    elif bound in ['DM','np']:
        mask = res[bound+'_trans'].mask
        mask_val = mask if np.prod(mask.shape)==1. else mask[p,y,n]

        return res[bound+'_trans'][p,y,n] if ~mask_val else None
    else:
        assert False, 'Not implemented!'
