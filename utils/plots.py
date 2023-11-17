import numpy as np

import matplotlib.pyplot as plt
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

from DM_theory.network import network

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
            ax.plot(res[x_key][:trans_idx],res['q'][0,a,:trans_idx],color=(col,0,0),label=r'$\displaystyle \alpha_0 = %g$'%res['alpha_0'][a])
        else:
            ax.plot(res[x_key],res['q'][0,a,],color=(col,0,0),label=r'$\displaystyle \alpha_0 = %g$'%res['alpha_0'][a])

    if approx:
        ## plotting approximations
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.1*steps)],res[x_key][:int(0.1*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (steps1 > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,ymax*1.1],
        # ylim=[0,20],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
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
        # ax.plot(res[x_key],res[x_key]**2,'k--',linewidth=0.5)
        # ax.plot(res[x_key][:int(0.012*steps)],res[x_key][:int(0.012*steps)]*res['rateWnt'][-1]/np.sqrt(2),'r--',linewidth=0.5)
        ax.plot(res[x_key],res['q_approx'][0,0,:],'k--')
        if (len(res[x_key]) > 1):
            ax.plot(res[x_key],res['q_approx'][0,2,:],'r--')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,5],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle q\,$[Hz$\displaystyle ^2$]'
    )

    if order:
        set_title(ax,order=order,title=r'the low $\displaystyle \bar{\nu}$ limit')


def plot_currents(ax,res,plt_para,bound='imp',idxs=None,approx=False,plot_options=['var','I'],order=1):
    
    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape
    print(steps1,steps)
            
    ## plot constant lines
    ax.axhline(0,c='k',ls='--',lw=0.5)

    if not idxs:
        idxs = range(steps1)

    if 'var' in plot_options:
    
        ## plot dark-matter transitions for one of the solutions
        trans_idx = get_trans_idx(res,'DM',2,0,0)
        if trans_idx:
            val = res[x_key][trans_idx]
            ax.axvline(val,color='k',ls=':')
            # ax.text(val+0.5,0.2,r'$\displaystyle \nu_{DM}$',fontsize=10)
        
        ## plotting temporal fluctuations
        trans_idx = get_trans_idx(res,bound,0,0,0)
        ax.plot(res[x_key][:trans_idx],res['sigma_V'][0,0,:trans_idx],color='k',ls='--',label=r'$\displaystyle \sigma_V$')
        ax.plot(res[x_key][trans_idx:],res['sigma_V'][0,0,trans_idx:],color='k',ls=':')
        
        ## plotting quenched fluctuations
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            trans_idx = get_trans_idx(res,bound,i,0,0)
            
            ax.plot(res[x_key][:trans_idx],res['alpha'][0,i,:trans_idx],color=(col,0,0),ls='-',label=r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$')
            ax.plot(res[x_key][trans_idx:],res['alpha'][0,i,trans_idx:],color=(col,0,0),ls=':')

        ax.text(res[x_key][int(steps/2)],res['sigma_V'][0,0,int(steps/2)]*1.1,r'$\displaystyle \sigma_{V_k}$',fontsize=10)
        ax.text(res[x_key][int(steps*2/3)],res['alpha'][0,0,int(steps*2/3)]*0.7,r'$\displaystyle \alpha_k = \sqrt{\alpha_{I_k}^2 + \alpha_0^2}$',fontsize=10)

        plt.setp(ax, ylim=[-0.02,0.15])
        
    if 'I' in plot_options:
        
        ## plotting currents
        for i,a in enumerate(idxs):
            col = i/float(len(idxs)-1)
            trans_idx = get_trans_idx(res,bound,i,0,0)
            
            ax.plot(res[x_key],-res['I_balance'][0,i,:],':',color=[0.7,0,0],lw=0.8)
            ax.plot(res[x_key][:trans_idx],-res['I_balance'][0,i,:trans_idx],'-',color=(col,0,0),label=r'$\displaystyle I_{balance}$',lw=0.8)

        ax.text(res[x_key][int(steps/2)],-res['I_balance'][0,0,int(steps/2)]*0.9,r'$\displaystyle \bar{I}_0-\Psi_0$',fontsize=10)

        plt.setp(ax,
            ylim=[-0.3,0.02]
        )

    if 'var' in plot_options and 'I' in plot_options:
        plt.setp(ax,
            ylim=[-0.3,0.3]
        )

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle I$'
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
        col = i/float(len(idxs)-1)
        # if a!=1:
        for p in range(Npop):
            trans_idx = get_trans_idx(res,bound,a,0,0)
            ax.plot(res[x_key],res['gamma'][p,a,:]**2,color=(col,0,0),ls=':')
            ax.plot(res[x_key][:trans_idx],res['gamma'][p,a,:trans_idx]**2,color=(col,0,0))

            if approx:
                ax.plot(res[x_key],res['gamma_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    if order:
        set_title(ax,order=order,title='$\displaystyle\quad$ DM exponent $\displaystyle \gamma$')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,7],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle \gamma^2$'
    )


def plot_chi(ax,res,plt_para,bound='imp',idxs=None,approx=False,order=1):

    x_key = plt_para['x']['key']
    Npop,steps1,steps = res['q'].shape

    if not idxs:
        idxs = range(steps1)

    for i,a in enumerate(idxs):
        col = i/float(len(idxs)-1)
        p=0
        # for p in range(Npop):
        trans_idx = get_trans_idx(res,bound,a,0,p)
        ax.plot(res[x_key],res['chi'][p,a,:],color=(col,0,0),ls=':')
        ax.plot(res[x_key][:trans_idx],res['chi'][p,a,:trans_idx],color=(col,0,0))

        if approx:
            ax.plot(res[x_key],res['chi_approx'][p,a,:]**2,color=(col,0,0),ls='--')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,15],
        xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',
        ylabel=r'$\displaystyle \chi$'
    )

    if order:
        set_title(ax,order=order,title='skewness coeff. $\displaystyle \chi$')


def plot_regions(ax,res,plt_para,order=1):

    # Npop,steps1,steps = res['q'].shape

    mask_DM = ~res['DM_trans'][0,:,0].mask
    mask_no_peak = ~res['np_trans'][0,:,0].mask
    mask_inc = ~res['inc_trans'][:,0].mask if np.any(res['inc_trans'][:,0].mask) else np.ones_like(res['inc_trans'][:,0],'bool')

    x_arr = res[plt_para['x']['key']]

    ax.plot(x_arr[res['DM_trans'][0,mask_DM,0]],res['alpha_0'][mask_DM],'r-',label=r'$\displaystyle \bar{\nu}_{DM}$')
    ax.plot(x_arr[res['np_trans'][0,mask_no_peak,0]],res['alpha_0'][mask_no_peak],'k--',label=r'$\displaystyle \bar{\nu}_{no\,peak}$')
    ax.plot(x_arr[res['inc_trans'][mask_inc,0]],res['alpha_0'][mask_inc],'k',label=r'$\displaystyle \bar{\nu}_{inc}$')
    #ax.plot(x_arr[res['imp_trans'][mask_imp,0]],res['alpha_0'][mask_imp],'k:',label=r'$\displaystyle \bar{\nu}_{imp}$')

    plt.setp(ax,
        xlim=[0,plt_para['x']['lim']],
        ylim=[0,0.12],
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


def plot_approx(ax_ex,ax_I,ax_alpha):
    # print()

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
    ax_ex.axhline(-np.sqrt(net.I_squared_nu(nu,nu**2,p)),color=[0.6,0.6,0.6],lw=2,ls='--')#,label=r"approx.: $\displaystyle q=\bar{\nu}^2$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_nu(nu,q,p)),'k-',label=r"solution for $\displaystyle \bar{\nu}$")
    ax_ex.plot(q,-np.sqrt(net.I_squared_q(nu,q,p)),'k--',label=r"solution for $\displaystyle q$")
    ax_ex.annotate(r'$\displaystyle q = \bar{\nu}^2$',xy=[nu**2,-np.sqrt(net.I_squared_nu(nu,nu**2,p))],xytext=[nu**2-22,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)
    ax_ex.annotate(r'$\displaystyle (\bar{\nu}^{\star},q^{\star})$',xy=[q[min_idx],-np.sqrt(net.I_squared_nu(nu,q[min_idx],p))],xytext=[nu**2+15,-0.18],arrowprops=dict(arrowstyle="->"),fontsize=10)

    plt.setp(ax_ex, ylim=ylims)
    ax_ex.legend(prop={'size':10},bbox_to_anchor=(1.2,0.0),loc='lower right',frameon=False)
    loc = ticker.MultipleLocator(base=0.1) # this locator puts ticks at regular intervals
    ax_ex.yaxis.set_major_locator(loc)
    plot_I(ax_I,res_approx,'k')
    plot_I(ax_I,res_approx_hetero,'r')

    plt.setp(ax_I,xlabel=r'$\displaystyle \bar{\nu}$ [Hz]',ylabel=r'$\displaystyle \frac{|\Delta I|}{\Psi_0 - I_0}$')

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


def get_displayString(key):
    
    if key == 'rateWnt':
        return r'$\displaystyle \bar{\nu}\,$[Hz]'
    elif key == 'alpha_0':
        return r'$\displaystyle \alpha_0$'
    elif key == 'tau_I':
        return r'$\displaystyle \tau_I\,$[ms]'
    elif key == 'eps':
        return r'$\displaystyle \varepsilon$'
    elif key == 'eta':
        return r'$\displaystyle \eta$'
    elif key == 'n':
        return r'$\displaystyle n$'
    elif key == 'Psi_0':
        return r'$\displaystyle \Psi_0$'
    else:
        return r''


def get_trans_idx(res,bound,y,n=0,p=0):

    if not ((bound+'_trans') in res.keys()):
        return None
    if bound in ['imp','inc']:
        mask = res[bound+'_trans'].mask
        mask_val = mask if np.prod(mask.shape)==1 else mask[y,n]
        
        mask = np.all(res[bound+'_trans'].mask)
        return res[bound+'_trans'][y,n] if ~mask_val else None
    elif bound in ['DM','np']:
        mask = res[bound+'_trans'].mask
        mask_val = mask if np.prod(mask.shape)==1. else mask[p,y,n]

        return res[bound+'_trans'][p,y,n] if ~mask_val else None
    else:
        assert False, 'Not implemented!'
