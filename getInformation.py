import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
import scipy.stats as stats
import matplotlib as mpl
from matplotlib.widgets import Slider

from netCDF4 import Dataset, stringtochar

from darkMatter import darkMatter
from plotting.statistics import *
from pythonCode.network import network

class DM_information:

    order=['rateWnt','alpha_0','tau_G','n','eta','eps','I_alpha','I_beta']
    rateWnt=1.
    alpha_0=0.01
    tau_G=0.005

    I_alpha=1.0
    I_beta=1.0

    arrays = {
        'rateWnt':[0,20],
        'alpha_0':[0,0.1],
        'tau_G':[0.005,0.1],
        'I_alpha':[0.,3.],
        'I_beta':[0.,3.],
    }

    options = {'mode_stats': 4}
    plt_para = {}

    def __init__(self,steps=100,
        save=0,file_format='png',
        rerun=False,compile=False):

        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')

        levs = 20
        self.plt_para['bnw'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=levs-1)
        self.plt_para['heat'] = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=levs-1)
        self.plt_para['bnw_regions'] = mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3)


        self.steps = steps
        self.set_arrays()

        self.net = network()
        # self.rate_arr, distr = self.net.distribution(1,1,steps=1000)

        self.rerun = rerun
        self.compile = compile

        ## plot preparation
        self.fig = plt.figure(figsize=(8,6),dpi=300)

        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())

        self.axes = {
            'I_fun': plt.axes([0.175,0.1,0.2,0.15]),
            'DM_phase': plt.axes([0.5,0.05,0.3,0.5],projection="3d"),
            'I_phase': plt.axes([0.5,0.55,0.3,0.5],projection="3d"),
            'distr': plt.axes([0.175,0.8,0.2,0.15]),
            'I_distr_cum': plt.axes([0.175,0.6,0.2,0.15])
        }
        self.axes['I_distr'] = self.axes['distr'].twinx()

        self.handles = {}

        self.plot_I()
        self.plot_DM()


        ## construct sliders
        slider_w = 0.02
        slider_h = 0.35
        axamp_alpha = plt.axes([0.02, .6, slider_w, slider_h])
        axamp_tau = plt.axes([0.06, .6, slider_w, slider_h])
        axamp_I_alpha = plt.axes([0.03, .1, slider_w, slider_h])
        axamp_I_beta = plt.axes([0.06, .1, slider_w, slider_h])
        axamp_rate = plt.axes([0.1, 0.6, slider_w, slider_h])

        self.slider = {}

        self.slider['alpha_0'] = Slider(axamp_alpha, r'$\displaystyle \alpha_0$', self.arrays['alpha_0'][0], self.arrays['alpha_0'][-1], valinit=self.alpha_0,orientation='vertical')
        self.slider['tau_G'] = Slider(axamp_tau, r'$\displaystyle \tau_G$', self.arrays['tau_G'][0], self.arrays['tau_G'][-1], valinit=self.tau_G,orientation='vertical')
        self.slider['I_alpha'] = Slider(axamp_I_alpha, r'$\displaystyle I_{\alpha}$', self.arrays['I_alpha'][0], self.arrays['I_alpha'][-1], valinit=self.I_alpha,orientation='vertical')
        self.slider['I_beta'] = Slider(axamp_I_beta, r'$\displaystyle I_{\beta}$', self.arrays['I_beta'][0], self.arrays['I_beta'][-1], valinit=self.I_beta,orientation='vertical')
        self.slider['rateWnt'] = Slider(axamp_rate, r'$\displaystyle \bar{\nu}$', self.arrays['rateWnt'][0], self.arrays['rateWnt'][-1], valinit=self.rateWnt,orientation='vertical')

        self.update_DM(0)
        self.update_I(0)

        # call update function on slider value change
        self.slider['alpha_0'].on_changed(self.update_DM)
        self.slider['tau_G'].on_changed(self.update_DM)
        self.slider['I_alpha'].on_changed(self.update_I)
        self.slider['I_beta'].on_changed(self.update_I)

        self.slider['rateWnt'].on_changed(self.update_all)


        plt.show()

        # ax_I_fun =
        # ax_phase =
        # ax_phase2 =
        # ax_2d_1 = plt.axes([0.75,0.6,0.2,0.15])
        # ax_2d_2 = plt.axes([0.75,0.4,0.2,0.15])
        #
        # ax_2d_3 = plt.axes([0.75,0.2,0.2,0.15])
        # ax_2d_4 = plt.axes([0.75,0.05,0.2,0.15])

    def set_arrays(self):

        self.arrays['alpha_0'] = np.linspace(self.arrays['alpha_0'][0],self.arrays['alpha_0'][-1],self.steps)
        self.arrays['tau_G'] = np.linspace(self.arrays['tau_G'][0],self.arrays['tau_G'][-1],self.steps)
        self.arrays['I_alpha'] = np.linspace(self.arrays['I_alpha'][0],self.arrays['I_alpha'][-1],self.steps)
        self.arrays['I_beta'] = np.linspace(self.arrays['I_beta'][0],self.arrays['I_beta'][-1],self.steps)
        self.arrays['rateWnt'] = np.linspace(self.arrays['rateWnt'][0],self.arrays['rateWnt'][-1],self.steps)

    def plot_DM(self):

        res = get_DM_phase(self.options,self.steps,rateWnt=self.rateWnt,I_alpha=self.I_alpha,I_beta=self.I_beta,rerun=self.rerun)

        x = np.linspace(stats.beta.ppf(0.001, self.I_alpha, self.I_beta),stats.beta.ppf(0.99, self.I_alpha, self.I_beta), 1000)
        y = stats.beta.pdf(x,self.I_alpha, self.I_beta)
        if ('I_fun' in self.handles):
            y = stats.beta.pdf(x,self.I_alpha,self.I_beta)
            self.handles['I_fun'].set_ydata(y)
        else:
            self.handles['I_fun'], = self.axes['I_fun'].plot(x,y)
            plt.setp(self.axes['I_fun'],
                xlabel=r'$\displaystyle \nu / \nu_{max}$',
                ylabel=r'$\displaystyle \beta(x,\alpha,\beta)$')

        plt.setp(self.axes['I_fun'],ylim=[0,max(1.5,np.nanpercentile(y,95)*1.1)])

        if ('DM_phase' in self.handles):
            try:
                self.handles['DM_phase'].remove()
            except:
                pass
        self.handles['DM_phase'] = plt_phase(self.axes['DM_phase'],res,self.plt_para,'alpha_0','tau_G')

        self.update_distr()

    def plot_I(self):

        res = get_I_phase(self.options,self.steps,rateWnt=self.rateWnt,alpha_0=self.alpha_0,tau_G=self.tau_G,rerun=self.rerun)

        idx_alpha_0 = self.get_idx('alpha_0')
        idx_tau_G = self.get_idx('tau_G')
        self.q = res['q'][0, idx_alpha_0, idx_tau_G]

        if ('I_phase' in self.handles):
            try:
                self.handles['I_phase'].remove()
            except:
                pass
        self.handles['I_phase'] = plt_phase(self.axes['I_phase'],res,self.plt_para,'I_alpha','I_beta')

        self.update_distr()


    def update_distr(self):
        idx_alpha_0 = self.get_idx('alpha_0')
        idx_tau_G = self.get_idx('tau_G')

        rate_arr, distr = self.net.distribution(self.rateWnt,self.q,steps=1000)

        if ('distr' in self.handles):
            self.handles['distr'].set_ydata(distr)
        else:
            self.handles['distr'], = self.axes['distr'].plot(rate_arr,distr,'k-')
        plt.setp(self.axes['distr'],xlim=[0.,0.2],ylim=[0,np.nanmax(distr)*1.1], ylabel=r'$\displaystyle \rho(\nu)$')

        y = rate_arr*distr*stats.beta.pdf(rate_arr/self.net.rate_max(),self.I_alpha,self.I_beta)
        if ('I_distr' in self.handles):
            self.handles['I_distr'].set_ydata(y)
        else:
            self.handles['I_distr'], = self.axes['I_distr'].plot(rate_arr,y,'r-')
        plt.setp(self.axes['I_distr'],xlim=[0.,0.5],ylim=[0,np.nanmax(y)*1.1], ylabel=r'$\displaystyle f(I(\nu))$')

        y = np.nancumsum(y)
        y /= y[-2]
        if ('I_distr_cum' in self.handles):
            self.handles['I_distr_cum'].set_ydata(y)
        else:
            self.handles['I_distr_cum'], = self.axes['I_distr_cum'].plot(rate_arr,np.nancumsum(y),'r-')
        plt.setp(self.axes['I_distr_cum'],xlim=[0.,0.5],ylim=[0,1],xlabel=r'$\displaystyle \nu / \nu_{max}$', ylabel=r'$\displaystyle F(I(\nu))$')


    def update_I(self,val):

        self.I_alpha = self.slider['I_alpha'].val
        self.I_beta = self.slider['I_beta'].val
        self.plot_DM()

        self.fig.canvas.draw_idle()


    def update_DM(self,val):

        self.alpha_0 = self.slider['alpha_0'].val
        self.tau_G = self.slider['tau_G'].val
        self.net.tau_I = self.tau_G
        self.net.alpha_0 = self.alpha_0

        self.plot_I()

        self.fig.canvas.draw_idle()

    def update_all(self,val):

        self.rateWnt = self.slider['rateWnt'].val
        self.update_DM(0)
        self.update_I(0)


    def get_idx(self,key):
        val = getattr(self, key)
        return np.nanargmin(abs(val-self.arrays[key]))



# def information(steps=100,
#     rateWnt=[0,20],alpha_0=[0,0.1],tau_G=[0.005,0.1],eps=[0.5],eta=[0.9],n=[0],I_alpha=[0.,2.5],I_beta=[0.,2.5],
#     order=['rateWnt','alpha_0','tau_G','n','eta','eps','I_alpha','I_beta'],
#     J=-1.,Npop=1,drive=0,
#     save=0,file_format='png',
#     rerun=False,compile=False):



def get_DM_phase(options,steps,I_alpha,I_beta,alpha_0=[0.,0.1],tau_G=[0.,0.1],rateWnt=[0.5],rerun=True,compile=False):

    options['order'] = ['alpha_0','tau_G','I_alpha','I_beta','rateWnt','n','eta','eps']
    options['rateWnt'] = [rateWnt]
    options['alpha_0'] = alpha_0
    options['tau_G'] = tau_G
    options['I_alpha'] = [I_alpha]
    options['I_beta'] = [I_beta]

    res = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)
    return res

def get_I_phase(options,steps,alpha_0,tau_G,rateWnt=[1.0],I_alpha=[0.,2.5],I_beta=[0.,2.5],rerun=True,compile=False):

    options['order'] = ['I_alpha','I_beta','alpha_0','tau_G','rateWnt','n','eta','eps']
    options['rateWnt'] = [rateWnt]
    options['alpha_0'] = [alpha_0]
    options['tau_G'] = [tau_G]
    options['I_alpha'] = I_alpha
    options['I_beta'] = I_beta

    return darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)



# def plt_I_fun(ax,alpha,beta):
#
#     x = np.linspace(stats.beta.ppf(0.01, alpha, beta),stats.beta.ppf(0.99, alpha, beta), 100)
#     y = stats.beta.pdf(x,alpha,beta)
#     h_p = ax.plot(x,y,label=r"$\displaystyle \alpha=%g$"%alpha)
#     plt.setp(ax,
#         xlabel=r'$\displaystyle \nu / \nu_{max}$',
#         ylabel=r'$\displaystyle \beta(x,\alpha,\beta)$',
#         ylim=[0,np.nanmedian(y)*2.1])
#
#     return h_p


def plt_phase(ax,res,plt_para,x_key,y_key):

    res['infoContent'][res['infoContent']==0] = np.nan
    idx_max1 = np.nanargmax(res['infoContent'][0,...],axis=0)
    idx_max2 = np.nanargmax(res['infoContent'][0,...],axis=1)

    col = get_col(res,plt_para)

    y_arr = res[y_key][:]
    x_arr = res[x_key][:]
    X,Y = np.meshgrid(x_arr,y_arr)
    h_p = ax.plot_surface(X,Y,res['infoContent'][0,...],facecolors=col,antialiased=False,shade=True)

    # ax.plot(res[options['order'][1]][idx_max1],res[options['order'][0]][range(steps)],res['infoContent'][0,range(steps),idx_max1]+0.1,'r',lw=3)
    # ax.plot(res[options['order'][1]][idx_max2],res[options['order'][0]][range(steps)],res['infoContent'][0,idx_max2,range(steps)]+0.001,'green',lw=3)

    plt.setp(ax,#zlim=[0,10],
        xlabel=r'$\displaystyle %s$'%x_key,
        ylabel=r'$\displaystyle %s$'%y_key,
        zlabel=r'$\displaystyle I_{\sum}$',
        zlim=[np.nanmin(res['infoContent']),np.nanmax(res['infoContent'])])


    # steps = len(x_arr)

    # h_p_2d_1 = []
    # h_p, = ax_2d[0].plot(x_arr,res['infoContent'][0,idx_max1,range(steps)],label=r'$\displaystyle I(%s^{max})$'%x_key)
    # h_p_2d_1.append(h_p)
    # h_p, = ax_2d[0].plot(y_arr,res['infoContent'][0,range(steps),idx_max2],label=r'$\displaystyle I(%s^{max})$'%y_key)
    # h_p_2d_1.append(h_p)
    # ax_2d[0].legend()
    # plt.setp(ax_2d[0],xlim=[0,0.1],
    #             xlabel='$\displaystyle %s \; or \; %s$'%(x_key,y_key))
    #
    #
    # h_p_2d_2 = []
    # h_p, = ax_2d[1].plot(x_arr,y_arr[idx_max1],label=r'$\displaystyle %s(%s^{max})$'%(y_key,x_key))
    # h_p_2d_2.append(h_p)
    # h_p, = ax_2d[1].plot(x_arr[idx_max2],y_arr,label=r'$\displaystyle %s(%s^{max})$'%(x_key,y_key))
    # h_p_2d_2.append(h_p)
    #
    # ax_2d[1].legend()
    # plt.setp(ax_2d[1],xlim=[0,0.1],
    #             xlabel='$\displaystyle %s$'%x_key)
    #
    # plt.suptitle(r'$\displaystyle \bar{\nu}=%.2fHz, \alpha = %.2f, \beta = %.2f$'%(res['rateWnt'][0],res['I_alpha'][0],res['I_beta'][0]))

    return h_p#, h_p_2d_1, h_p_2d_2


def get_col(res,plt_para):

    mask_dark_matter = (res['gamma'][0,...]**2 < 1)

    normalize = mcolors.Normalize(vmin=0, vmax=3)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['heat'])
    col_chi = s_map.to_rgba(res['chi'][0,...])

    normalize = mcolors.Normalize(vmin=0, vmax=2)
    s_map = cm.ScalarMappable(norm=normalize, cmap=plt_para['bnw'])
    col_gamma = s_map.to_rgba(res['gamma'][0,...])

    col = col_chi
    col[mask_dark_matter,:] = col_gamma[mask_dark_matter,:]

    col[0,0,:] =[np.nan,np.nan,np.nan,np.nan]

    return col
