import numpy as np
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
# from scipy.io import netcdf
import os, math
#import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
#from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
#from matplotlib._png import read_png
#from matplotlib.cbook import get_sample_data
import matplotlib.image as mpimg
from scipy.interpolate import griddata

from netCDF4 import Dataset


def sharkfins(steps=10,rateWnt=None,alpha_0=None,tau_G=None,n=None,eps=None,eta=None,Npop=1,drive=0,tau_M=10.,tau_A=5.,kappa=1,mode_calc="exact",mode_stats=0,plot_ax3D=False,save=0,file_format='png',compile=True,rerun=False):

## stats:
####    0: sharkfins
####    1: rate_dependence stats
####    2: borders phase space
####    3: compare exact vs. approx (single)
####    4: KL-phase-space (costly!)

    plt.rc('text', usetex=True)
    plt.rc('font', family='sans-serif')
    #mpl.rcParams['axes.labelsize'] = 'large'
    plt.rcParams['xtick.labelsize'] = 12
    #plt.rcParams['xtick.major.size'] = 2
    plt.rcParams['ytick.labelsize'] = 12
    #plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12        #### how to get proper and same fontsizes with no specified xticklabels and with specified ones?
    #steps = steps + 1
    if (drive == 1):
        assert Npop>1, "for drive == 1, you need at least two populations"

    max_steps = 0

    inpara = {}
    #### alpha vs nu (sharkfins)
    if tau_G:
    #assert len(tau_G) == len(inpara['rateWnt']), "tau_G should have the same length as rateWnt"
        inpara['tau_G'] = tau_G#[list(np.array(entry)*tau_M/0.01) for entry in tau_G] #tau_G
        max_steps = max(len(inpara['tau_G']),max_steps)
    else:
        #inpara['tau_G'] = [[0.005],[0.02],[0.04]]
        inpara['tau_G'] = [[5],[10],[20],[40],[60]]
    # [print(tau_list) for tau_list in inpara['tau_G']]
    # inpara['tau_G'] = [tau/1000 for tau in tau_list for tau_list in inpara['tau_G']]
    # print('tauG: ',inpara['tau_G'])

    if rateWnt:
        inpara['rateWnt'] = rateWnt#[list(np.array(entry)*0.01/tau_M) for entry in rateWnt] #rateWnt
        max_steps = max(len(inpara['rateWnt']),max_steps)
    else:
        #inpara['rateWnt'] = [[0,20],[0,10],[0,7]]
        inpara['rateWnt'] = [[0,rate[0]] for rate in (1./(2.*math.pi*np.sqrt(tau_M/1000.*np.array(inpara['tau_G'])/1000.)))]#[[0,20],[0,15],[0,10],[0,7],[0,6]]
        # print("rates: ", inpara['rateWnt'])

    if n:
        #assert len(n) == len(inpara['rateWnt']), "n should have the same length as rateWnt"
        inpara['n'] = n
        max_steps = max(len(inpara['n']),max_steps)
    else:
        #inpara['n'] = [[0],[0],[0]]
        inpara['n'] = [[0],[0],[0],[0],[0]]

    if alpha_0:
        #assert len(alpha_0) == len(inpara['rateWnt']), "alpha_0 should have the same length as rateWnt"
        inpara['alpha_0'] = alpha_0
        max_steps = max(len(inpara['alpha_0']),max_steps)
    else:
        #inpara['alpha_0'] = [[0,0.16],[0,0.08],[0,0.05]]
        inpara['alpha_0'] = [[0,0.16],[0,0.12],[0,0.08],[0,0.05],[0,0.035]]

    if eps:
        #assert len(eps) == len(inpara['rateWnt']), "eps should have the same length as rateWnt"
        inpara['eps'] = eps
        max_steps = max(len(inpara['eps']),max_steps)
    else:
        if Npop>1:
            inpara['eps'] = [[0.],[0.1],[0.3],[0.5],[math.sqrt(0.5)]]
        else:
            inpara['eps'] = [[0.],[0.],[0.],[0.],[0.]]

    if eta:
        #assert len(eta) == len(inpara['rateWnt']), "eta should have the same length as rateWnt"
        inpara['eta'] = eta
        max_steps = max(len(inpara['eta']),max_steps)
    else:
        if Npop>1:
            inpara['eta'] = [[0.],[0.2],[0.5],[0.8],[0.9]]
        else:
            inpara['eta'] = [[0.],[0.],[0.],[0.],[0.]]

    for para in ['rateWnt','n','alpha_0','tau_G','eps','eta']:
        if (max_steps > len(inpara[para])):
            inpara[para]*=max_steps

    #### tau vs nu
    # inpara['n'] = [[0],[0],[0]]
    # inpara['alpha_0']= [[0.01],[0.02],[0.04]]
    # inpara['tau_G'] = [[0,100],[0,100],[0,100]]
    # inpara['rateWnt'] = [[0,20],[0,20],[0,20]]

    #### alpha vs tau
    # inpara['n'] = [[0],[0],[0]]
    # inpara['alpha_0'] = [[0,0.1],[0,0.1],[0,0.1]]
    # inpara['tau_G'] = [[0,100],[0,100],[0,100]]
    # inpara['rateWnt'] = [[1],[2],[5]]

    #results = {}

    levs = range(20)
    plot_para = {
        'multi': len(inpara['rateWnt']) > 1,
        'cb_plotted': False,
        'bnw': mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=len(levs)-1),
        'heat': mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=len(levs)-1),
        'bnw_regions': mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3),
    }

    if (not mode_stats):
        if plot_para['multi']:
            if (Npop == 1):
                fig, ax = plt.subplots(len(inpara['rateWnt'])//2+1,2,figsize=(7.5,1+2.3*(len(inpara['rateWnt'])//2+1)))
            # if (Npop == 2):
                # fig, ax = plt.subplots(2,len(inpara['rateWnt']),figsize=(7.5,1+2.3*(len(inpara['rateWnt'])//2+1)))

        else:
            if (Npop == 1):
                fig = plt.figure(figsize=(6,5))
            if (Npop == 2):
                fig = plt.figure(figsize=(6,10))

    if (compile):
        compile_code()


    # set model parameters
    model = set_model(Npop,tau_A,tau_M,kappa,drive)

    fileModel = './data/modPara.nc'
    fileSim = './data/simPara.nc'
    fileResults = './data/results.nc'

    results = []

    for i in range(len(inpara['rateWnt'])):

        if plot_para['multi']:
            ax_now = ax[i//2,i%2]
        else:
            if (Npop == 1):
                ax_now = plt.axes([0.12,0.1,0.7,0.8])
            if (Npop == 2):
                ax_now = plt.axes([0.12,0.55-p*0.45,0.7,0.4])

        # if (rateWnt == None):
            # print("tau_G: ", inpara['tau_G'][i][0])
            # print("max rate: ", 1./(2.*math.pi*np.sqrt(model['tau_M']/1000.*inpara['tau_G'][i][0]/1000.)))

        para_list = ['alpha_0','n','tau_G','rateWnt','eps','eta']
        simulation, ax_list, ax_labels, const_labels, sv_str = set_simulation(inpara,i,steps,para_list,mode_calc,mode_stats)

        fileResults = './data/results_%dsteps_const_%s=%g_mode=%d_drive=%d.nc' % (steps,para_list[ax_list[2]],simulation[para_list[ax_list[2]]],mode_stats,drive)
        if not os.path.exists(fileResults) or rerun:
            write_input(model,fileModel,simulation,fileSim)
            run_str = './sharkfins ' + fileModel + ' ' + fileSim + ' ' + fileResults

            # print("Running: ",run_str)
            os.system(run_str)


        # obtain data from simulation
        ncid = Dataset(fileResults, "r")

        results.append({})
        print('reading results from %s ...'%fileResults)
        for key in ncid.variables:
            # print(key)
            # print(ncid.variables[key])
            results[i][key] = ncid.variables[key][:]
        ncid.close()

        if ((Npop==2) and not (mode_stats==1)):
            return ncid, para_list, inpara, ax_list, ax_labels, simulation

        if (steps == 1):
            return results

            plt.figure()
            plt.plot(p_range,cdf_theory)
            plt.show(block=False)

        else:

            if ((Npop==1) and (mode_stats==0)):
                for p in range(Npop):

                    plot_fins(ax_now,simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],results[i]['gamma'][p,...],results[i]['chi'][p,...],results[i]['regions'][p,...],plot_para)
                    set_axes(ax_now,simulation,para_list[ax_list[1]],para_list[ax_list[0]])

                    if (not plot_para['multi']):
                        ax_now.set_title('%s' % (const_labels[0]),fontsize=12)

    # plt.show(block=False)
    # return results

    if (steps > 1):
        if (mode_stats == 0):
            if (len(inpara['rateWnt']) > 1):

                big_ax = plt.axes([0.1,0.1,0.8,0.85])
                big_ax.set_facecolor('none')
                big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
                big_ax.spines['top'].set_visible(False)
                big_ax.spines['right'].set_visible(False)
                big_ax.spines['bottom'].set_visible(False)
                big_ax.spines['left'].set_visible(False)
                big_ax.set_xlabel(ax_labels[1],fontsize=12)
                big_ax.set_ylabel(ax_labels[0],fontsize=12)

            if (plot_para['multi']):
                plot_3d(fig,ax[-1,-1],inpara,results,para_list[ax_list[1]],para_list[ax_list[2]],para_list[ax_list[0]],ax_labels,const_labels,from_disk=not plot_ax3D)

            if (Npop==1):


                #ax[i/2,i%2].set_xlabel(ax_labels[1],fontsize=12)
                        #ax[i/2,i%2].set_ylabel(ax_labels[0],fontsize=12)

                if not plot_para['multi']:# and (steps>1)):
                    ax.set_xlabel(ax_labels[1],fontsize=16)
                    ax.set_ylabel(ax_labels[0],fontsize=16)

                #plt.suptitle('%s, %s' % (const_labels[0],const_labels[0]),fontsize=12)
                #plt.rcParams['xtick.labelsize'] = 16
                #mpl.rcParams['xtick.labelsize'] = 16
                if save:
                    sv_name = '../paper draft/inhib only/pics/transp_hai_steps=%d_%s_%s_drive%d.%s' % (steps,para_list[ax_list[0]],para_list[ax_list[1]],drive,file_format)
                    plt.savefig(sv_name,dpi=300)
                    print('Figure saved as "%s"' % sv_name)

                plt.show(block=False)

    elif (mode_stats in [1,3]):
        results['rate_max'] = (2*math.pi*np.sqrt(simulation['tau_G'][0]*model['tau_M']))**(-1)

        stps_alpha = len(results['gamma'])
        results['nu_c_idx'] = np.zeros((stps_alpha,2))
        results['nu_c'] = np.zeros(stps_alpha)

        for a in range(stps_alpha):

            if any(results['gamma'][a] > 1):
                #print results['gamma'][a]
                softDM = np.where(results['gamma'][a]>1)[0]
                #print softDM
                if len(softDM):
                    results['nu_c_idx'][a][0] = softDM[0]
                    results['nu_c_idx'][a][1] = softDM[-1]
                    results['nu_c'][a] = results['rateWnt'][0,results['nu_c_idx'][a][0]]
                else:
                    results['nu_c_idx'][a][0] = 0
                    results['nu_c_idx'][a][1] = -1
                    results['nu_c'][a] = np.nan
                #print results['nu_c_idx'][a]
                #results['nu_c'][a] = results['rateWnt'][0,results['nu_c_idx'][a][0]]
                #ax[1,0].plot([nu_c[a],nu_c[a]],[-1,results['gamma'][a,nu_c_idx]**2-1],'k--',linewidth=0.5)
            else:
                results['nu_c_idx'][a] = 0
                results['nu_c'][a] = np.nan
                #results['nu_c_idx'][a][1] = 0
        # print(results['nu_c'])
        return(results)

    elif (mode_stats == 2):
        results['rateWnt'] = simulation['rateWnt']
        results['alpha_0'] = simulation['alpha_0']
        return results

    elif (mode_stats == 3):
        return results

    elif (mode_stats == 4):
        return results,simulation

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

def compile_code():
    ## prepare code by compiling, if needed
    modes = '-g -Wall -ansi'

    libs = ' -lgsl -lgslcblas -lnetcdf -std=c++11'
    os.system('g++ %s -o sharkfins sharkfins.cpp %s' % (modes,libs))

    ## ---------------------------------------
    ## this is left just as an example of how to set paths (as I always forget...)
    #codepath = 'g++ -g -Wall -ansi  -o sharkfins -I/home/wollex/libs/gsl/include -I/home/wollex/libs/netcdf/include  sharkfins.cpp -L/home/wollex/libs/gsl/lib -lgsl -lgslcblas -L/home/wollex/libs/netcdf/lib -lnetcdf_c++ -lnetcdf -std=c++11'

    #libdir = ['gsl','netcdf']
    #pathlib = '/home/wollex/libs/'

    #buildL = ''
    #buildI = ''
    #for lib in libdir:
    #print lib
    #buildL += '-L' + pathlib + lib + '/lib '
    #buildI += '-I' + pathlib + lib + '/include '
    ## ---------------------------------------


def plot_fins(ax,x_arr,y_arr,gamma,chi,regions,plot_para):

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
    plot_para['bnw'].set_bad('k',0.)
    plot_para['heat'].set_bad('k',0.)
    pgamma = ax.pcolormesh(x_arr,y_arr,plot_gamma,cmap=plot_para['bnw'],vmin=0,vmax=2)
    pchi = ax.pcolormesh(x_arr,y_arr,plot_chi,cmap=plot_para['heat'],vmin=0,vmax=3)
    pregions = ax.pcolormesh(x_arr,y_arr,plot_regions,cmap=plot_para['bnw_regions'],vmin=2,vmax=3)
    pimplausible = ax.pcolormesh(x_arr,y_arr,plot_implausible,cmap=plot_para['bnw_regions'],vmin=1,vmax=3,alpha=0.4)

    plot_colorbar(pchi,pgamma,plot_para)

    if not plot_para['multi']:

        #ax[i/2,i%2].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_no_peak,cmap=bnw,vmin=-1,vmax=2)
        #ax[i/2,i%2].pcolormesh(simulation[para_list[ax_list[1]]],simulation[para_list[ax_list[0]]],plot_inconsistent,cmap=bnw,vmin=-1,vmax=2)

        mask_DM = ~np.isnan(results['DM_trans'])
        mask_no_peak = ~np.isnan(results['no_peak_trans'])
        mask_inc = ~np.isnan(results['inc_trans'])

        ax.plot(results['DM_trans'][mask_DM],y_arr[mask_DM],'r-',linewidth=2,label=r'$\displaystyle \bar{\nu}_{DM}$')
        ax.plot(results['no_peak_trans'][mask_no_peak],y_arr[mask_no_peak],'k--',linewidth=2,label=r'$\displaystyle \bar{\nu}_{no\,peak}$')
        ax.plot(results['inc_trans'][mask_inc],y_arr[mask_inc],'k',linewidth=2,label=r'$\displaystyle \bar{\nu}_{inc}$')
        ax.plot(results['nu_implausible'][mask_inc],y_arr[mask_inc],'k:',linewidth=2,label=r'$\displaystyle \bar{\nu}_{imp}$')

        # ax.set_xlabel(ax_labels[1],fontsize=12)
        # ax.set_ylabel(ax_labels[0],fontsize=12)

    ax.set_xlim([0,x_arr[-1]])
    ax.set_ylim([0,y_arr[-1]])


def set_axes(ax,simulation,x_key,y_key):

    for key in [x_key,y_key]:
        ticks = 11
        if key == 'n':
            array = np.linspace(0,1,6)
        elif key == 'alpha_0':
            if simulation[key][-1] > 0.1:
                ticks = 6
            array = np.linspace(0,0.2,ticks)
        elif key == 'tau_G':
            if simulation[key][-1] > 50:
                ticks = 6
            array = np.linspace(0,100,ticks)
        elif key == 'rateWnt':
            if simulation[key][-1] > 10:
                ticks = 6
            array = np.linspace(0,20,ticks)

        if key==x_key:
            try:
                x_array = array[:np.where(array>simulation[key][-1])[0][0]]
            except:
                x_array = array[:]
            # ax.set_xticks(x_array)
            # ax.set_x
        else:
            try:
                y_array = array[:np.where(array>simulation[key][-1])[0][0]]
            except:
                y_array = array[:]
            # ax.set_yticks(y_array)


def set_model(Npop,tau_A,tau_M,kappa,drive):

    # defining model constants (hard coded & accessed by programm call)
    model = {}
    model['Npop'] = Npop
    model['tau_A'] = tau_A/1000.  				# the excitatory time constant for AMPA currents
    model['tau_N'] = 200/1000.				    # the excitatory time constant for NMDA currents
    model['tau_M'] = tau_M/1000				    # the membrane time constant of all neurons
    #model['tau_G'] = inpara['tau_G']           # the inhibitory time constant for GABA currents
    #!!! nochmal checken: tau-G

    model['kappa'] = kappa
    model['drive'] = drive

    if (drive > 0):
        model['tau_0'] = 1.
        model['J_0'] = 1.
        model['K_0'] = 1
    else:
        model['tau_0'] = 0.
        model['J_0'] = 0.
        model['K_0'] = 0

    return model


def set_simulation(inpara,i,steps,para_list,mode_calc,mode_stats):

    simulation = {}
    for key in inpara:
        assert type(inpara[key][i]) == list, 'Please specify all iteration parameters as lists!'
        if len(inpara[key][i]) == 2:
            simulation[key] = np.linspace(inpara[key][i][0],inpara[key][i][-1],steps+1)[1:]#.astype('float')
        elif len(inpara[key][i]) > 2:
            simulation[key] = np.array(inpara[key][i])
        else:
            simulation[key] = np.array(inpara[key][i]).astype('float')

    simulation['tau_G'] /= 1000

    ax_list = []
    ax_labels = []
    const_labels = []
    sv_str = 'steps=%d_' % steps

    idx = 0

    sv_str += '_iter'
    for key in para_list:
        if len(simulation[key]) > 1:
            sv_str += '_%s=%g' % (key,simulation[key][-1])
            ax_list.append(idx)

            if key=='eps':
                ax_labels.append(r'$\displaystyle \varepsilon$')
            elif key == 'eta':
                ax_labels.append(r'$\displaystyle \eta$')
            elif key == 'n':
                ax_labels.append(r'$\displaystyle b$')
            elif key == 'alpha_0':
                ax_labels.append(r'$\displaystyle \alpha_0$')
            elif key == 'tau_G':
                ax_labels.append(r'$\displaystyle \tau_G\,$[ms]')
            elif key == 'rateWnt':
                ax_labels.append(r'$\displaystyle \bar{\nu}\,$[Hz]')

        idx += 1

    idx = 0
    sv_str += 'const'
    for key in para_list:

        if len(simulation[key]) == 1:
            if (len(inpara[key]) == 1):
                idx += 1
                continue
            elif (inpara[key][0] == inpara[key][1]):
                idx += 1
                continue

            ax_list.append(idx)
            sv_str += '_%s=%g' % (key,simulation[key][0])
            if key == 'eps':
                const_labels.append(r'$\displaystyle \varepsilon = %g$' % simulation[key])
            elif key == 'eta':
                const_labels.append(r'$\displaystyle \eta = %g$' % simulation[key])
            elif key == 'n':
                const_labels.append(r'$\displaystyle n = %g$' % simulation[key])
            elif key == 'alpha_0':
                const_labels.append(r'$\displaystyle \alpha_0 = %g$' % simulation[key])
            elif key == 'tau_G':
                const_labels.append(r'$\displaystyle \tau_I = %g\,$ms' % simulation[key])
            elif key == 'rateWnt':
                const_labels.append(r'$\displaystyle \bar{\nu} = %g\,$Hz' % simulation[key])

        idx += 1

    simulation['mode_calc'] = int(mode_calc=='approx')
    simulation['mode_stats'] = int(mode_stats)

    return simulation, ax_list, ax_labels, const_labels, sv_str

def write_input(model,fileModel,simulation,fileSim):

    #write model parameters to netcdf file
    ncid = Dataset(fileModel, "w");#, format="NETCDF4")
    ncid.createDimension('one',1)
    for key in model.keys():
        if type(model[key])==int:
            varType = 'i'   #np.dtype('int32').char
        else:
            varType = 'f8'  #np.dtype('float64').char

        Var = ncid.createVariable(key,varType,('one',))
        Var[:] = model[key]
    ncid.close()


    #write simulation parameters to netcdf file
    ncid = Dataset(fileSim, "w")
    ncid.createDimension('one',1)
    for key in simulation:
        if type(simulation[key]) == np.ndarray:
            varDim = key + 'Sz'
            ncid.createDimension(varDim,len(simulation[key]))
            varType = simulation[key][0].dtype.char
        else:
            varDim = 'one'
            if type(simulation[key])==int:
                varType = 'i'   #np.dtype('int32').char
            else:
                varType = 'f8'  #np.dtype('float64').char
            # var_type = type(simulation[key])
            # var_type = 'f8'

        Var = ncid.createVariable(key,varType,(varDim,))
        Var[:] = simulation[key]

    # Var = ncid.createVariable('mode_calc','i',('one',))
    # if mode_calc == 'exact':
    #     Var[:] = 0
    # elif mode_calc == 'approx':
    # Var[:] = mode_calc=='approx'

    # Var = ncid.createVariable('mode_stats','i',('one',))
    # Var[:] = mode_stats

    ncid.close()

def plot_3d(fig,ax,inpara,results,px,py,pz,ax_labels,const_labels,from_disk=False):

    print("plotting 3D axis")

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

        DM_bound = np.zeros((len(inpara['rateWnt']),steps)) * np.nan
        DM_bound_max = np.zeros((len(inpara['rateWnt']),2)) * np.nan
        imp_bound = np.zeros((len(inpara['rateWnt']),steps)) * np.nan
        inc_bound = np.zeros((len(inpara['rateWnt']),steps)) * np.nan

        data_3D_inc = np.zeros([len(inpara['rateWnt']),steps,3]) * np.nan
        data_3D_DM = np.zeros([len(inpara['rateWnt']),steps,3]) * np.nan

        fig.delaxes(ax)
        ax = fig.add_subplot(len(inpara['rateWnt'])//2+1,2,len(inpara['rateWnt'])+1,projection='3d')

        # fig2 = plt.figure()
        # ax3D = fig2.add_subplot(111,projection='3d')

        print('axes: ',px,py,pz)
        cross_point = np.zeros(len(results),'int');

        for i,res in enumerate(results):

            print('px: ',px)
            print('py: ',py)
            print()
            print('pz: ',pz)

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
                cross_point[i] = np.where(data_3D_DM[i,:,0]>data_3D_inc[i,:,0])[0][0]
            else:
                cross_point[i] = np.where(~np.isnan(data_3D_DM[i,:,0]))[0][-1]

        i_arr = np.arange(len(results))
        ax.plot(data_3D_inc[i_arr,cross_point,0],data_3D_inc[i_arr,cross_point,1],data_3D_inc[i_arr,cross_point,2],'k-',lw=2)
        ax.plot(data_3D_inc[i_arr,0,0],data_3D_inc[i_arr,0,1],data_3D_inc[i_arr,0,2],'k-',lw=2)
        ax.plot(data_3D_DM[i_arr,0,0],data_3D_DM[i_arr,0,1],data_3D_DM[i_arr,0,2],'k-',lw=2)

        for i in range(len(results)-1):

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
        ax.dist = 8

        ax.set_title('%s' % (const_labels[0]),fontsize=12)

        # ax.set_xticks(np.linspace(inpara[px][0][0],inpara[px][0][1],5))
        # ax.tick_params(axis='x', pad=-5)
        # ax.set_yticks(inpara[py][0],inpara[py][-1])#np.linspace(inpara[para_list[ax_list[2]]][0][0],inpara[para_list[ax_list[2]]][-1][0],5))
        # ax.set_yticklabels(['%d'%(entry[0]) for entry in inpara[pz]],fontsize=12)#np.linspace(inpara[para_list[ax_list[2]]][0][0],inpara[para_list[ax_list[2]]][-1][0],5))
        # ax.tick_params(axis='y', pad=-3)
        # ax.set_zticks(np.linspace(inpara[pz][0][0],inpara[pz][0][1],5))

        # ax.set_xlim(x_plane)
        # ax.set_ylim(y_plane)

        # ax.set_zlim([0,0.12])
        ax.set_xlabel(ax_labels[1],labelpad=-5)
        ax.set_ylabel(const_labels[0],labelpad=-5)#const_labels[0])
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel(ax_labels[0],rotation=90)


def plot_colorbar(pchi,pgamma,plot_para):

    if not plot_para['cb_plotted']:
        plot_para['cb_plotted'] ^= True;

        axcb1 = plt.axes([0.85,0.35,0.03,0.6])
        axcb2 = plt.axes([0.85,0.1,0.03,0.2])

        axcb1.tick_params(axis='both', which='major', labelsize=12)
        axcb2.tick_params(axis='both', which='major', labelsize=12)

        plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,7))
        plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,10),ticks=np.linspace(0,1,3))

        #axcb2.set_yticks(np.linspace(1,0,3))
        #axcb2.set_yticklabels(np.linspace(1,0,3))
        axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=12)
        axcb2.set_ylabel(r'$\displaystyle \gamma^2$',fontsize=12)

        plt.subplots_adjust(left=0.12, bottom=0.1, right=0.8, top=0.95, wspace=0.35, hspace=0.3)
