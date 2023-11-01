import numpy as np
import matplotlib.pyplot as plt
import math

from darkMatter import darkMatter

from general.plot_statistics import *
from general.utils import set_plot_params

# def sharkfins(steps=10,rateWnt=None,alpha_0=None,tau_G=None,n=None,eps=None,eta=None,Npop=1,drive=0,tau_M=10.,tau_A=5.,kappa=1,mode_calc="exact",mode_stats=0,plot_ax3D=False,save=0,file_format='png',compile=True,rerun=False):
def sharkfins(steps=10,Npop=1,plot_ax3D=True,save=0,file_format='png',rerun=False,compile=False):

## stats:
####    0: sharkfins
####    1: rate_dependence stats
#######    2: borders phase space
####    2: compare exact vs. approx (single)
####    3: KL-phase-space (costly!)

    steps = steps + 1       # correct for removal of first item

    ## general plot setup
    set_plot_params()

    inpara = {}

    #### alpha vs nu (sharkfins)
    para_order = ['rateWnt','alpha_0','tau_G','n','eps','eta']
    tau_M = 0.01
    inpara['tau_G'] = [[0.005],[0.01],[0.02],[0.04],[0.1]]
    inpara['alpha_0'] = [[0,0.16],[0,0.12],[0,0.08],[0,0.05],[0,0.035]]
    inpara['rateWnt'] = [[0,rate[0]] for rate in (1./(2.*math.pi*np.sqrt(tau_M*np.array(inpara['tau_G']))))]#[[0,20],[0,15],[0,10],[0,7],[0,6]]

    if Npop==1:
        inpara['eta'] = [[0.],[0.],[0.],[0.],[0.]]
        inpara['eps'] = [[0.],[0.],[0.],[0.],[0.]]
        inpara['n'] = [[0],[0],[0],[0],[0]]
    else:
        inpara['eta'] = [[0.],[0.2],[0.5],[0.8],[0.9]]
        inpara['eps'] = [[0.],[0.1],[0.3],[0.5],[math.sqrt(0.5)]]
        inpara['n'] = [[0],[0],[0],[0],[0]]


    #### tau vs nu
    # para_order = ['rateWnt','tau_G','alpha_0','n','eps','eta']
    # inpara['n'] = [[0],[0],[0]]
    # inpara['tau_G'] = [[0,100],[0,100],[0,100]]
    # inpara['alpha_0']= [[0.01],[0.02],[0.04]]
    # inpara['rateWnt'] = [[0,20],[0,20],[0,20]]

    #### alpha vs tau
    # para_order = ['tau_G','alpha_0','rateWnt','n','eps','eta']
    # inpara['n'] = [[0],[0],[0]]
    # inpara['tau_G'] = [[0,100],[0,100],[0,100]]
    # inpara['alpha_0'] = [[0,0.1],[0,0.1],[0,0.1]]
    # inpara['rateWnt'] = [[1],[2],[5]]

    ## brief sanity checks on provided parameters
    sim_steps = 0
    for key in para_order:
        assert type(inpara[key]) == list, 'Please specify all parameters as lists of lists!'
        assert type(inpara[key][0]) == list, 'Please specify all parameters as lists of lists!'
        sim_steps = max(sim_steps,len(inpara[key]))

    for key in ['rateWnt','n','alpha_0','tau_G','eps','eta']:
        if len(inpara[key]) == 1:
            inpara[key]*=sim_steps
        else:
            assert sim_steps == len(inpara[key]), 'Make sure all parameters are specified as having the same length (%s has length %d / %d)' % (key,len(inpara[key]),sim_steps)


    levs = range(20)
    plt_para = {
        'multi': sim_steps > 1,
        'ax_label': [],
        'const_label': []
    }

    if plt_para['multi']:
        fig, ax = plt.subplots(sim_steps//2+1,
                    2,
                    figsize=(7.5,1+2.3*(sim_steps//2+1)))
        big_ax = fig.add_axes([0.1,0.1,0.8,0.85])
        big_ax.set_facecolor('none')
        big_ax.tick_params(labelcolor='none',top=False,bottom=False,left=False,right=False)
        big_ax.spines['top'].set_visible(False)
        big_ax.spines['right'].set_visible(False)
        big_ax.spines['bottom'].set_visible(False)
        big_ax.spines['left'].set_visible(False)
    else:
        fig = plt.figure(figsize=(6,5) if Npop==1 else (6,10))


    # return
    results = []
    for i in range(sim_steps):
        # print(i)
        options = {}
        for key in para_order:
            options[key] = inpara[key][i]

        results.append(darkMatter(steps=steps,options=options,rerun=rerun,compile=compile))

        set_labels(results[i],para_order,plt_para)

        for p in range(Npop):
            if plt_para['multi']:
                ax_now = ax[i//2,i%2]
            else:
                if (Npop == 1):
                    ax_now = plt.axes([0.12,0.1,0.7,0.8])
                if (Npop == 2):
                    ax_now = plt.axes([0.12,0.55-p*0.45,0.7,0.4])

            for p in range(Npop):
                plot_fins(ax_now,results[i][para_order[0]],results[i][para_order[1]],results[i]['gamma'][p,...],results[i]['chi'][p,...],results[i]['regions'][p,...],plt_para)
                set_axes(ax_now,results[i],para_order[0],para_order[1])

                if (not plt_para['multi']):
                    ax_now.set_title('%s' % (plt_para['const_label'][0]),fontsize=12)

    if plt_para['multi']:
        big_ax.set_xlabel(plt_para['ax_label'][1],fontsize=12)
        big_ax.set_ylabel(plt_para['ax_label'][0],fontsize=12)
    else:
        ax_now.set_xlabel(plt_para['ax_label'][1],fontsize=12)
        ax_now.set_ylabel(plt_para['ax_label'][0],fontsize=12)

    # if (steps > 1):
        # if (mode_stats == 0):
    if (plt_para['multi']):
        plot_3d(fig,ax[-1,-1],results,para_order,plt_para,from_disk=not plot_ax3D)

    plt.suptitle('%s' % (plt_para['const_label'][0]),fontsize=12)
    if save:
        sv_name = '../figures/shark_steps=%d_%s_%s_drive%d.%s' % (steps,para_order[0],para_order[1],drive,file_format)
        plt.savefig(sv_name,dpi=300)
        print('Figure saved as "%s"' % sv_name)

    plt.show(block=False)

    # if (mode_stats in [1,3]):
    #     results['rate_max'] = (2*math.pi*np.sqrt(simulation['tau_G'][0]*model['tau_M']))**(-1)
    #
    #     stps_alpha = len(results['gamma'])
    #     results['nu_c_idx'] = np.zeros((stps_alpha,2))
    #     results['nu_c'] = np.zeros(stps_alpha)
    #
    #     for a in range(stps_alpha):
    #
    #         if any(results['gamma'][a] > 1):
    #             #print results['gamma'][a]
    #             softDM = np.where(results['gamma'][a]>1)[0]
    #             #print softDM
    #             if len(softDM):
    #                 results['nu_c_idx'][a][0] = softDM[0]
    #                 results['nu_c_idx'][a][1] = softDM[-1]
    #                 results['nu_c'][a] = results['rateWnt'][0,results['nu_c_idx'][a][0]]
    #             else:
    #                 results['nu_c_idx'][a][0] = 0
    #                 results['nu_c_idx'][a][1] = -1
    #                 results['nu_c'][a] = np.nan
    #             #print results['nu_c_idx'][a]
    #             #results['nu_c'][a] = results['rateWnt'][0,results['nu_c_idx'][a][0]]
    #             #ax[1,0].plot([nu_c[a],nu_c[a]],[-1,results['gamma'][a,nu_c_idx]**2-1],'k--',linewidth=0.5)
    #         else:
    #             results['nu_c_idx'][a] = 0
    #             results['nu_c'][a] = np.nan
    #             #results['nu_c_idx'][a][1] = 0
    #     # print(results['nu_c'])
    #     return(results)
    #
    # elif (mode_stats == 2):
    #     results['rateWnt'] = simulation['rateWnt']
    #     results['alpha_0'] = simulation['alpha_0']
    #     return results
    #
    # elif (mode_stats == 3):
    #     return results
    #
    # elif (mode_stats == 4):
    #     return results,simulation

    return results
    #print "test4"

  #os.remove(fileResults)

def do_all():
  haifische(steps=500,rate=[0,20],alpha_0=[0,0.16],tau_G=[0.005],n=[0])
  haifische(steps=500,rate=[0,15],alpha_0=[0,0.12],tau_G=[0.01],n=[0])
  haifische(steps=500,rate=[0,10],alpha_0=[0,0.08],tau_G=[0.02],n=[0])
  haifische(steps=500,rate=[0,7],alpha_0=[0,0.05],tau_G=[0.04],n=[0])
  haifische(steps=500,rate=[0,6.5],alpha_0=[0,0.035],tau_G=[0.06],n=[0])

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
        pos2 = [pos1.x0 - 0.175, pos1.y0-0.025,  pos1.width + 0.2, pos1.height + 0.075]

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
        pos2 = [pos1.x0 - 0.1, pos1.y0-0.05,  pos1.width + 0.2, pos1.height + 0.075]

        fig.delaxes(ax)
        ax = plt.axes(pos2,projection='3d')

        cross_point = np.zeros(len(results),'int');

        maxPara = {'x':0,'y':0,'z':0}
        for i,res in enumerate(results):

            maxPara['x'] = max(maxPara['x'],np.nanmax(res[px]))
            maxPara['y'] = max(maxPara['y'],np.nanmax(res[py]))
            maxPara['z'] = max(maxPara['z'],np.nanmax(res[pz]))

            inc_bound[i,:] = res['inc_trans'][0,:]
            DM_bound[i,:] = res['DM_trans'][0,:]

            data_3D_DM[i,:,0] = res['DM_trans']
            data_3D_DM[i,:,1] = res[py]
            data_3D_DM[i,:,2] = res[pz]

            data_3D_inc[i,:,0] = res['inc_trans']
            data_3D_inc[i,:,1] = res[py]
            data_3D_inc[i,:,2] = res[pz]

            ax.plot(data_3D_DM[i,:,0],data_3D_DM[i,:,1],data_3D_DM[i,:,2],'k-',lw=2)
            ax.plot(data_3D_inc[i,:,0],data_3D_inc[i,:,1],data_3D_inc[i,:,2],'k-',lw=2)

            # print('crossing point: @')
            cross = np.where(data_3D_DM[i,:,0]>data_3D_inc[i,:,0])[0]
            if len(cross):
                cross_point[i] = np.where((data_3D_DM[i,:,0]>data_3D_inc[i,:,0]) & (~np.isnan(data_3D_DM[i,:,0]) | np.isnan(data_3D_inc[i,:,0])))[0][0]
            else:
                cross_point[i] = np.where(~np.isnan(data_3D_DM[i,:,0]))[0][-1]

        print(maxPara)
        i_arr = np.arange(sim_steps)
        ax.plot(data_3D_inc[i_arr,cross_point,0],data_3D_inc[i_arr,cross_point,1],data_3D_inc[i_arr,cross_point,2],'k-',lw=2)
        ax.plot(data_3D_inc[i_arr,0,0],data_3D_inc[i_arr,0,1],data_3D_inc[i_arr,0,2],'k-',lw=2)
        ax.plot(data_3D_DM[i_arr,0,0],data_3D_DM[i_arr,0,1],data_3D_DM[i_arr,0,2],'k-',lw=2)

        for i in range(sim_steps-1):

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
        ax.azim = -100.
        ax.elev = 20
        # ax.dist = 8

        # ax.set_xticks(np.linspace(inpara[px][0][0],inpara[px][0][1],5))
        # ax.tick_params(axis='x', pad=-5)
        ax.set_yticks(np.linspace(0,maxPara['y'],3))#np.linspace(inpara[para_order[ax_list[2]]][0][0],inpara[para_order[ax_list[2]]][-1][0],5))
        # ax.set_yticklabels(['%d'%(entry[0]) for entry in inpara[pz]],fontsize=12)#np.linspace(inpara[para_order[ax_list[2]]][0][0],inpara[para_order[ax_list[2]]][-1][0],5))
        ax.set_yticklabels(np.linspace(0,maxPara['y'],3), ha='right', rotation_mode="anchor",rotation=20)
        ax.tick_params(axis='y', pad=-2)
        # ax.set_zticks(np.linspace(inpara[pz][0][0],inpara[pz][0][1],5))

        ax.set_xlim([0,maxPara['x']*1.0])
        ax.set_ylim([0,maxPara['y']*1.2])
        ax.set_zlim([0,maxPara['z']*1.1])

        ax.set_xlabel(plt_para['ax_label'][0],labelpad=5)
        ax.set_ylabel(plt_para['ax_label'][2],labelpad=8)#const_labels[0])
        ax.set_zlabel(plt_para['ax_label'][1],ha='left',rotation_mode='anchor',rotation=90,labelpad=5)
        ax.zaxis.set_rotate_label(False)
