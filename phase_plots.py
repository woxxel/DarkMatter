import numpy as np
import matplotlib.pyplot as plt
import math

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *

# def sharkfins(steps=10,rateWnt=None,alpha_0=None,tau_G=None,n=None,eps=None,eta=None,Npop=1,drive=0,tau_M=10.,tau_A=5.,kappa=1,mode_calc="exact",mode_stats=0,plot_ax3D=False,save=0,file_format='png',compile=True,rerun=False):
def phase_plots(steps=100,Npop=1,plot_ax3D=False,save=0,file_format='png',rerun=False,compile=False):

## stats:
####    0: sharkfins
####    1: rate_dependence stats
#######    2: borders phase space
####    2: compare exact vs. approx (single)
####    3: KL-phase-space (costly!)

    steps = steps + 1       # correct for removal of first item

    ## general plot setup
    # plot parameters
    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    # plt.rcParams['font.family'] = ['Tahoma','Verdana']
    plt.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    title_offset = -0.0
    fig,ax = plt.subplots(2,2,figsize=(7.5,7.5),dpi=300)
    plt_para = {
        'ax_label': [],
        'const_label': []
    }

    p = 0
    options = {
        'order': ['rateWnt','alpha_0','tau_G','n','eta','eps'],
        'rateWnt': [0,20],
        'alpha_0': [0,0.2],
        'mode_stats': 0
    }
    results_1 = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    plot_fins(ax[0,0],results_1[options['order'][0]],results_1[options['order'][1]],results_1['gamma'][p,...],results_1['chi'][p,...],results_1['regions'][p,...],plt_para)
    plt.setp(ax[0,0],ylabel=r'$\displaystyle \alpha_0$',xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]')
    ax[0,0].tick_params(labelbottom=False, labeltop=True, labelleft=True, labelright=False,
                    bottom=True, top=True, left=True, right=True)
    ax[0,0].xaxis.set_label_position('top')
    set_title(ax[0,0],1,r'$\displaystyle \tau_G = 5\,$ms',(title_offset,1.15),10)

    options = {
        'order': ['tau_G','alpha_0','rateWnt','n','eta','eps'],
        'tau_G': [0,0.1],
        'alpha_0': [0,0.2],
        'rateWnt': [1.],
        'mode_stats': 0
    }
    results_2 = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    plot_fins(ax[0,1],1000*results_2[options['order'][0]],results_2[options['order'][1]],results_2['gamma'][p,...],results_2['chi'][p,...],results_2['regions'][p,...],plt_para)
    plt.setp(ax[0,1],xlabel=r'$\displaystyle \tau_G\,$[ms]',ylabel=r'$\displaystyle \alpha_0$')
    ax[0,1].tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=True,
    bottom=True, top=True, left=True, right=True)
    ax[0,1].xaxis.set_label_position('top')
    ax[0,1].yaxis.set_label_position('right')
    set_title(ax[0,1],2,r'$\displaystyle \bar{\nu} = %g\,$Hz'%options['rateWnt'][0],(title_offset,1.15),10)

    options = {
        'order': ['rateWnt','tau_G','alpha_0','n','eta','eps'],
        'rateWnt': [0,20],
        'tau_G': [0,0.1],
        'alpha_0': [0.04],
        'mode_stats': 0
    }
    results_3 = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    plot_fins(ax[1,0],results_3[options['order'][0]],1000*results_3[options['order'][1]],results_3['gamma'][p,...],results_3['chi'][p,...],results_3['regions'][p,...],plt_para)
    plt.setp(ax[1,0],xlabel=r'$\displaystyle \bar{\nu}\,$[Hz]',ylabel=r'$\displaystyle \tau_G\,$[ms]')
    ax[1,0].tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                    bottom=True, top=True, left=True, right=True)
    set_title(ax[1,0],3,r'$\displaystyle \alpha_0 = %g\,$'%options['alpha_0'][0],(title_offset,1.2),20)


    #### alpha vs nu (sharkfins)
    inpara = {}
    tau_M = 0.01
    inpara['order'] = ['rateWnt','alpha_0','tau_G','n','eps','eta']
    inpara['tau_G'] = [[0.005],[0.007],[0.01],[0.015],[0.02],[0.03],[0.04],[0.05],[0.06]]#,[0.08],[0.1]]
    inpara['alpha_0'] = [[0,0.16],[0,0.16],[0,0.12],[0,0.12],[0,0.08],[0,0.08],[0,0.05],[0,0.05],[0,0.035]]#,[0,0.035],[0,0.035]]
    inpara['rateWnt'] = [[0,rate[0]] for rate in (1./(2.*math.pi*np.sqrt(tau_M*np.array(inpara['tau_G']))))]#[[0,20],[0,15],[0,10],[0,7],[0,6]]
    inpara['n'] = [[0]]
    inpara['eps'] = [[0]]
    inpara['eta'] = [[0]]


    ## brief sanity checks on provided parameters
    sim_steps = 0
    for key in inpara['order']:
        assert type(inpara[key]) == list, 'Please specify all parameters as lists of lists!'
        assert type(inpara[key][0]) == list, 'Please specify all parameters as lists of lists!'
        sim_steps = max(sim_steps,len(inpara[key]))

    for key in inpara['order']:
        if len(inpara[key]) == 1:
            inpara[key]*=sim_steps
        else:
            assert sim_steps == len(inpara[key]), 'Make sure all parameters are specified as having the same length (%s has length %d / %d)' % (key,len(inpara[key]),sim_steps)

    results_3D = []
    for i in range(sim_steps):
        # print(i)
        options = {}
        for key in inpara['order']:
            options[key] = inpara[key][i]

        results_3D.append(darkMatter(steps,options,rerun=rerun,compile=compile))
    # return results_3D
    set_labels(results_3D[0],inpara['order'],plt_para)

    ax[-1,-1] = plot_3d(fig,ax[-1,-1],results_3D,inpara['order'],plt_para,from_disk=False)
    # set_title(ax[-1,-1],3,r'phase space',(-0.1,0.3),-50)

    xx, zz = np.meshgrid(results_1['rateWnt'], results_1['alpha_0'])
    yy = xx*0 + results_1['tau_G'][0]*1000
    ax[-1,-1].plot_surface(xx, yy, zz, alpha=0.2, color='r')
    mask = ~np.isnan(results_1['DM_trans'][0,:])
    results_1['DM_trans'][0,results_1['DM_trans'][0,:]==0] = np.NaN
    results_1['inc_trans'][0,results_1['inc_trans'][0,:]==0] = np.NaN
    ax[-1,-1].plot(results_1['DM_trans'][0,:],np.ones(steps)*5,results_1['alpha_0'],'r-',lw=2)
    ax[-1,-1].plot(results_1['inc_trans'][0,mask],np.ones(mask.sum())*5,results_1['alpha_0'][mask],'r-',lw=2)

    yy, zz = np.meshgrid(results_2['tau_G']*1000, results_2['alpha_0'])
    xx = yy*0 + results_2['rateWnt'][0]
    results_2['DM_trans'][0,results_2['DM_trans'][0,:]==0] = np.NaN
    ax[-1,-1].plot_surface(xx, yy, zz, alpha=0.2, color='g')
    ax[-1,-1].plot(np.ones(steps)*2,1000*results_2['DM_trans'][0,:],results_2['alpha_0'],'g-',lw=2)
    # ax[-1,-1].plot(np.ones(steps)*2,1000*results_2['inc_trans'][0,:],results_2['alpha_0'],'g-',lw=1)

    xx, yy = np.meshgrid(results_3['rateWnt'], results_3['tau_G']*1000)
    zz = xx*0 + results_3['alpha_0'][0]
    ax[-1,-1].plot_surface(xx, yy, zz, alpha=0.2, color='b')
    mask = ~np.isnan(results_3['DM_trans'][0,:])
    results_3['DM_trans'][0,results_3['DM_trans'][0,:]==0] = np.NaN
    results_3['inc_trans'][0,results_3['inc_trans'][0,:]==0] = np.NaN
    ax[-1,-1].plot(results_3['DM_trans'][0,:],1000*results_3['tau_G'],np.ones(steps)*0.04,'b-',lw=2)
    ax[-1,-1].plot(results_3['inc_trans'][0,mask],1000*results_3['tau_G'][mask],np.ones(mask.sum())*0.04,'b-',lw=2)

    fig.tight_layout(pad=3.0)
    if save:
        sv_name = './figures/phase_plots.%s' % (file_format)
        plt.savefig(sv_name,format=file_format,dpi=600)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)
    return results_2

def cm2inch(value):
  return value/2.54


def set_axes(ax,results,x_key,y_key):

    for key in [x_key,y_key]:
        ticks = 11
        if key == 'n':
            array = np.linspace(0,1,6)
        elif key == 'alpha_0':
            if results[key][-1] > 0.1:
                ticks = 6
            array = np.linspace(0,0.2,ticks)
        elif key == 'tau_G':
            if results[key][-1] > 50:
                ticks = 6
            array = np.linspace(0,100,ticks)
        elif key == 'rateWnt':
            if results[key][-1] > 10:
                ticks = 6
            array = np.linspace(0,20,ticks)

        if key==x_key:
            try:
                x_array = array[:np.where(array>results[key][-1])[0][0]]
            except:
                x_array = array[:]
            # ax.set_xticks(x_array)
            # ax.set_x
        else:
            try:
                y_array = array[:np.where(array>results[key][-1])[0][0]]
            except:
                y_array = array[:]
            # ax.set_yticks(y_array)



def set_labels(results,para_order,plt_para):

    for key in para_order:

        if key=='eps':
            plt_para['ax_label'].append(r'$\displaystyle \varepsilon$')
        elif key == 'eta':
            plt_para['ax_label'].append(r'$\displaystyle \eta$')
        elif key == 'n':
            plt_para['ax_label'].append(r'$\displaystyle b$')
        elif key == 'alpha_0':
            plt_para['ax_label'].append(r'$\displaystyle \alpha_0$')
        elif key == 'tau_G':
            plt_para['ax_label'].append(r'$\displaystyle \tau_G\,[ms]$')
        elif key == 'rateWnt':
            plt_para['ax_label'].append(r'$\displaystyle \bar{\nu}\,[Hz]$')

        if len(results[key]) == 1:
            plt_para['const_label'].append(plt_para['ax_label'][-1][:-1] + ' = %g$'%results[key][0])
            break


def plot_3d(fig,ax,results,para_order,plt_para,from_disk=False):

    print("plotting 3D axis")

    px,pz,py = para_order[:3]
    if from_disk:
        pos1 = ax.get_position() # get the original position
        pos2 = [pos1.x0 - 0.175, pos1.y0-0.025,  pos1.width + 0.25, pos1.height + 0.175]

        fig.delaxes(ax)
        ax = plt.axes(pos2)

        if (pz == 'tau_G'):
            img=mpimg.imread("./figures/diagram_tau_3D_adjusted.png")
            print("tau")
        if (pz == 'rateWnt'):
            img=mpimg.imread("./figures/diagram_nu_3D_adjusted.png")
            print("rate")
        if (pz == 'alpha_0'):
            img=mpimg.imread("./figures/diagram_alpha_3D_adjusted.png")
            print("alpha")

        ax.imshow(img)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.setp(ax,xticks=[],yticks=[])

    else:
        Npop,steps,_ = results[0]['gamma'].shape
        sim_steps = len(results)
        DM_bound = np.zeros((sim_steps,steps)) * np.nan
        DM_bound_max = np.zeros((sim_steps,2)) * np.nan
        imp_bound = np.zeros((sim_steps,steps)) * np.nan
        inc_bound = np.zeros((sim_steps,steps)) * np.nan

        data_3D_inc = np.zeros([sim_steps,steps,3]) * np.nan
        data_3D_DM = np.zeros([sim_steps,steps,3]) * np.nan

        pos1 = ax.get_position() # get the original position
        pos2 = [pos1.x0 - 0.1, pos1.y0-0.05,  pos1.width + 0.25, pos1.height+0.125]

        fig.delaxes(ax)
        ax = plt.axes(pos2,projection='3d')
        ax.set_zorder(-10)

        cross_point = np.zeros(len(results),'int');

        maxPara = {'x':0,'y':0,'z':0}
        for i,res in enumerate(results):

            maxPara['x'] = max(maxPara['x'],np.nanmax(res[px]))
            maxPara['y'] = max(maxPara['y'],np.nanmax(res[py]))
            maxPara['z'] = max(maxPara['z'],np.nanmax(res[pz]))

            inc_bound[i,:] = res['inc_trans'][0,:]
            DM_bound[i,:] = res['DM_trans'][0,:]

            data_3D_DM[i,:,0] = res['DM_trans']
            data_3D_DM[i,:,1] = 1000*res[py]
            data_3D_DM[i,:,2] = res[pz]

            data_3D_inc[i,:,0] = res['inc_trans']
            data_3D_inc[i,:,1] = 1000*res[py]
            data_3D_inc[i,:,2] = res[pz]

            if not i%2:
                ax.plot(data_3D_DM[i,:,0],data_3D_DM[i,:,1],data_3D_DM[i,:,2],'k-',lw=1)
                ax.plot(data_3D_inc[i,:,0],data_3D_inc[i,:,1],data_3D_inc[i,:,2],'k-',lw=1)

            # print('crossing point: @')
            cross = np.where(data_3D_DM[i,:,0]>data_3D_inc[i,:,0])[0]
            if len(cross):
                cross_point[i] = np.where((data_3D_DM[i,:,0]>data_3D_inc[i,:,0]) & (~np.isnan(data_3D_DM[i,:,0]) | np.isnan(data_3D_inc[i,:,0])))[0][0]
            else:
                cross_point[i] = np.where(~np.isnan(data_3D_DM[i,:,0]))[0][-1]

            data_3D_inc[i,cross_point[i]+1:,0] = np.NaN
        # print(data_3D_DM[:,:,2])
        # print(data_3D_inc[:,:,2])
        # print(data_3D_DM[:,:,2] > )
        # print(maxPara)

        i_arr = np.arange(0,sim_steps,1)

        ax.plot(data_3D_inc[i_arr,cross_point[i_arr],0],data_3D_inc[i_arr,cross_point[i_arr],1],data_3D_inc[i_arr,cross_point[i_arr],2],'k-',lw=1)
        ax.plot(data_3D_inc[i_arr,0,0],data_3D_inc[i_arr,0,1],data_3D_inc[i_arr,0,2],'k-',lw=1)
        ax.plot(data_3D_DM[i_arr,0,0],data_3D_DM[i_arr,0,1],data_3D_DM[i_arr,0,2],'k-',lw=1)

        for i in range(sim_steps-1):
            # print(i)

            X = np.concatenate((data_3D_DM[i,:,0],np.flipud(data_3D_DM[i+1,:,0])))
            mask = ~np.isnan(X)
            Y = np.concatenate((data_3D_DM[i,:,1],np.flipud(data_3D_DM[i+1,:,1])))
            Z = np.concatenate((data_3D_DM[i,:,2],np.flipud(data_3D_DM[i+1,:,2])))
            X = np.pad(X[mask],1,mode='wrap')
            Y = np.pad(Y[mask],1,mode='wrap')
            Z = np.pad(Z[mask],1,mode='constant',constant_values=0)
            ax.plot_trisurf(X,Y,Z,color='grey',alpha=0.5)

            X = np.concatenate((data_3D_inc[i,:,0],np.flipud(data_3D_inc[i+1,:,0])))
            mask = ~np.isnan(X)
            Y = np.concatenate((data_3D_inc[i,:,1],np.flipud(data_3D_inc[i+1,:,1])))
            Z = np.concatenate((data_3D_inc[i,:,2],np.flipud(data_3D_inc[i+1,:,2])))
            X = np.pad(X[mask],1,mode='wrap')
            Y = np.pad(Y[mask],1,mode='wrap')
            Z = np.pad(Z[mask],1,mode='constant',constant_values=0)
            ax.plot_trisurf(X,Y,Z,color=(0.5,0.5,0.5,0.5),antialiased=True,linewidth=0,edgecolor='none',shade=True)
        ax.azim = -130.
        ax.elev = 15
        # ax.dist = 8

        # ax.set_xticks(np.linspace(inpara[px][0][0],inpara[px][0][1],5))
        # ax.tick_params(axis='x', pad=-5)
        # ax.set_yticklabels(['%d'%(entry[0]) for entry in inpara[pz]],fontsize=12)#np.linspace(inpara[para_order[ax_list[2]]][0][0],inpara[para_order[ax_list[2]]][-1][0],5))
        ax.set_yticks(np.linspace(0,1000*maxPara['y'],3))#np.linspace(inpara[para_order[ax_list[2]]][0][0],inpara[para_order[ax_list[2]]][-1][0],5))
        ax.set_yticklabels(['%d'%i for i in np.linspace(0,1000*maxPara['y'],3)], ha='right', rotation_mode="anchor",rotation=20)
        ax.set_zticks(np.linspace(0,0.1,3))#np.linspace(inpara[para_order[ax_list[2]]][0][0],inpara[para_order[ax_list[2]]][-1][0],5))
        # ax.set_zticklabels(['%d'%i for i in np.linspace(0,1000*maxPara['y'],3)], ha='right', rotation_mode="anchor",rotation=20)
        ax.tick_params(axis='y', pad=-2)
        # ax.set_zticks(np.linspace(inpara[pz][0][0],inpara[pz][0][1],5))

        ax.set_xlim([0,maxPara['x']*1.0])
        ax.set_ylim([0,1000*maxPara['y']*1.2])
        # ax.set_zlim([0,maxPara['z']*1.1])
        ax.set_zlim([0,0.14])

        ax.set_xlabel(plt_para['ax_label'][0],labelpad=5)
        ax.set_ylabel(plt_para['ax_label'][2],labelpad=8)#const_labels[0])
        ax.set_zlabel(plt_para['ax_label'][1],ha='left',rotation_mode='anchor',rotation=90,labelpad=5)
        ax.zaxis.set_rotate_label(False)
    return ax
