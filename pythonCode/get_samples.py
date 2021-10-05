import numpy as np
#from numpy.ma import masked_array
#import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os,math, imp
import matplotlib as mpl

#imp.load_source('read_mat', 'read_mat.py')
#from read_mat import *
import matplotlib.gridspec as gridspec


def plot_poisson(lam):

    plt.figure()
    plt.plot(np.arange(lam*5),get_poisson(lam)*(lam*5))
    plt.plot(np.arange(2*lam*5)/2.,get_poisson(lam*2)*(lam*2*5))
    plt.plot(np.arange(5*lam*5)/5.,get_poisson(lam*5)*(lam*5*5))

    plt.show(block=False)

    print(np.sum(get_poisson(lam)))
    print(np.sum(get_poisson(lam*2)))
    print(np.sum(get_poisson(lam*5)))

def get_poisson(lam):

    k_max = lam*5
    p = np.zeros(k_max)
    for i in range(k_max):
        p[i] = poisson_distr(lam,i)

    return p


def poisson_distr(lam,k):   # only single values allowed

    p = np.exp(k*np.log(lam)-lam-math.lgamma(k+1)) # numerically stable computation (avoiding overflow from power and factorial)

    return p


#def analyze_measure(fileName,prior='plain'):


    #data = mat2ncdf(fileName)


    ## specify the computation parameters
    #computation = {}

    #computation['p_theory'] = 0
    #computation['p_theory_hist'] = 0
    #if computation['p_theory_hist']:
        #assert computation['p_theory'], "To generate a discretized theoretical firing rate distribution, enable computation of theoretical computation!"
    #computation['draw_from_theory'] = 0                         # number of different draws from the theoretical distribution
    #computation['draw_finite_time'] = 0                         # number of finite time measurements per draw from theoretical distribution
    #computation['process_data'] = 1

    #computation['prior'] = prior


    ## specify file names
    #fileMeasures = './data/mesData.nc'
    #fileComp = './data/comPara.nc'
    #fileResults = './data/results.nc'

    #print "write netcdf to: ", fileMeasures

    #ncid = netcdf.netcdf_file(fileMeasures,'w')

    #ncid.createDimension('one',1)
    #ncid.createDimension('NSz',data['N'])

    #Var_T = ncid.createVariable('T','f',('one',))
    #Var_T[:] = data['T']

    #Var_N_AP = ncid.createVariable('N_AP','d',('NSz',))
    #Var_N_AP[:] = data['N_AP']

    #Var_rate = ncid.createVariable('rates','f',('NSz',))
    #Var_rate[:] = data['rates']

    #ncid.close()


    #ncid = netcdf.netcdf_file(fileComp,'w')
    #ncid.createDimension('one',1)

    #for key in computation.keys():
        #if key == 'prior':
            #ncid.createDimension('priorSz',len(computation['prior']))
            #Var = ncid.createVariable(key,'c',('priorSz',))
        #else:
            #Var = ncid.createVariable(key,'d',('one',))
        #Var[:] = computation[key]
    #ncid.close()

    #run_str = './single_distr ' + '1 ' + fileMeasures + ' ' + fileComp + ' ' + fileResults
    #os.system(run_str)


    #ncid = netcdf.netcdf_file(fileResults,'r',mmap=None)

    ##print ncid.variables.keys()

    #data['p_range'] = ncid.variables['p_range'][:]
    #data['p_bayes_est'] = ncid.variables['p_bayes_est'][:]
    ##data['rate_T'] = ncid.variables['rate_T'][:]

    ##plt.figure()
    #plt.plot(data['p_range'],data['p_bayes_est'])
    #plt.xlim([0,0.1])
    #plt.ylim([0,150])
    #plt.xlabel('firing rate (Hz)')
    #plt.ylabel('p(nu)')

    #plt.show(block=False)

    #return data



def get_samples_from_theory(rate=1,alpha_0=0,tau_G=0.005,tau_M=0.010,n=0.,eps=0,eta=0.9,Npop=1,drive=0,N=100,T=120,n_bin=50,prior="plain",p_theory=1,dftheory=1,dftime=1,mode_calc="exact",mode_stats=1,save=0,plot=0):

    simulation = {}

    simulation['n'] = n
    simulation['alpha_0'] = alpha_0
    simulation['tau_G'] = tau_G
    simulation['rateWnt'] = rate
    simulation['eps'] = eps
    simulation['eta'] = eta



    # defining model constants (hard coded & accessed by programm call)
    model = {}

    model['Npop'] = Npop

    model['tau_A'] = 0.005				    # the excitatory time constant for AMPA currents
    model['tau_N'] = 0.200				    # the excitatory time constant for NMDA currents
    model['tau_M'] = tau_M				    # the membrane time constant of all neurons

    #model['tau_G'] = tau_G				    # the inhibitory time constant for GABA currents
    #model['tau_G'] = tau_G
    #model['tau_G'] = 0.010/tau_M*tau_G
    #print "tau_G = ", model['tau_G']
    #model['n'] = n

    model['kappa'] = 1
    #model['eta'] = 0.9
    #model['eps'] = eps

    #model['rateWnt'] = rate
    #model['alpha_0'] = alpha_0

    model['drive'] = drive

    if (drive > 0):
        model['tau_0'] = 1
        model['J_0'] = 1
        model['K_0'] = 1
    else:
        model['tau_0'] = 0
        model['J_0'] = 0
        model['K_0'] = 0


    # specify the computation parameters
    computation = {}

    computation['n_bin'] = n_bin                             # number of bins for histogram construction
    computation['border'] = 1/(2*math.pi*np.sqrt(simulation['tau_G']*model['tau_M']))                         # upper border of AP numbers for calculations (hard coded, but should always be where p(nu)->0)
    computation['T'] = T                                     # duration of measurement

    computation['p_theory'] = p_theory                      # 0 = dont, 1 = calculate gamma, delta & chi, 2 = calculate full pdf (if possible)
    computation['p_theory_hist'] = 0
    if computation['p_theory_hist']:
        assert computation['p_theory'], "To generate a discretized theoretical firing rate distribution, enable computation of theoretical computation!"
    computation['draw_from_theory'] = dftheory                      # number of different draws from the theoretical distribution
    computation['draw_finite_time'] = dftime                     # number of finite time measurements per draw from theoretical distribution
    computation['process_data'] = 0

    if computation['draw_from_theory']:
        computation['N'] = N                                 # number of neurons measured
        computation['seed_theory'] = np.random.randint(1,2**52-1,dftheory)   # seed for random number generator to draw measurements
        computation['prior'] = prior


    if computation['draw_finite_time']:
        assert computation['draw_from_theory'], "To draw finite time samples, draw from theory has to be enabled!"
        computation['seed_time'] = np.random.randint(1,2**52-1,(dftheory*dftime))   # seed for random number generator to draw measurements


    ### define code paths
    if not os.path.exists('./data'):
        os.mkdir('./data')

    # specify file names
    fileModel = './data/modPara.nc'
    fileSim = './data/simPara.nc'
    fileComp = './data/comPara.nc'
    fileResults = './data/results.nc'

    #if not os.path.exists(fileResults):
    #print 'saving model parameters in %s ...' % fileNet
    ncid = netcdf.netcdf_file(fileModel,'w')
    ncid.createDimension('one',1)
    for key in model.keys():
        Var = ncid.createVariable(key,np.dtype('float32').char,('one',))
        Var[:] = model[key]
    ncid.close()

    #print 'saving computation parameters in %s ...' % fileSim
    ncid = netcdf.netcdf_file(fileComp,'w')
    ncid.createDimension('one',1)
    ncid.createDimension('two',2)

    ncid.createDimension('sims',dftheory*dftime)

    for key in computation.keys():
        if key == 'prior':
            ncid.createDimension('priorSz',len(computation['prior']))
            Var = ncid.createVariable(key,'c',('priorSz',))
        elif 'seed' in key:
            dim_seed = len(computation[key])
            if not (str(dim_seed) in ncid.dimensions.keys()):
                ncid.createDimension(str(dim_seed),dim_seed)
            Var = ncid.createVariable(key,'d',(str(dim_seed),))
        else:
            Var = ncid.createVariable(key,'d',('one',))
        Var[:] = computation[key]
    ncid.close()


    ncid = netcdf.netcdf_file(fileSim,'w')
    ncid.createDimension('one',1)

    for key in simulation:
        var_dim = key + 'Sz'
        ncid.createDimension(var_dim,1)
        var_type = type(simulation[key])

        if not type(simulation[key]) == np.ndarray:
            var_dim = 'one'
            var_type = type(simulation[key])
            var_type = 'd'
        Var = ncid.createVariable(key,var_type,(var_dim,))
        Var[:] = simulation[key]

    Var = ncid.createVariable('mode_calc','i',('one',))
    if mode_calc == 'exact':
        Var[:] = 0
    elif mode_calc == 'approx':
        Var[:] = 1

    Var = ncid.createVariable('mode_stats','i',('one',))
    Var[:] = mode_stats

    ncid.close()

    print(simulation)




    ### compile code if not yet done (convert to doing this via make file, which is executed here)
    if not os.path.exists('single_distr'):

        # define compiling variables
        mode = '-g -Wall -ansi'
        output = '-o single_distr'
        include = '-I /usr/local/include'
        libs = '-L /usr/local/lib -lgsl -lgslcblas -lnetcdf_c++ -lnetcdf'

        os.system('g++ %s %s %s single_distr.cpp %s' % (mode,output,include,libs))

    ### run code to generate distribution

    if computation['draw_from_theory']:
        print("Obtaining %d*%d samples from N=%d neurons " % (computation['draw_from_theory'],computation['draw_finite_time'],computation['N']))
    if computation['draw_finite_time']:
        print("with T=%d" % computation['T'])
    run_str = './single_distr ' + '0 ' + fileModel + ' ' + fileSim + ' ' + fileComp + ' ' + fileResults
    #print run_str
    os.system(run_str)


    ### read results from distribution

    x_lim = rate*5.
    l_H = x_lim/(n_bin)
    d_nu = 1./T

    results = {}
    #results['hist_range'] = np.linspace(0,x_lim,n_bin+1)

    #print 'reading from file %s ...' % fileResults
    ncid = netcdf.netcdf_file(fileResults,'r',mmap=None)

    if computation['draw_from_theory']:
        results['rate_inf'] = ncid.variables['rate_inf'][:]    # rate randomly generated
    if computation['draw_finite_time']:
        results['rate_T'] = ncid.variables['rate_T'][:]        # rate during finite time (from poisson drawn n_AP with rate_inf in time T)

    ###### is it possible to obtain stem-function (cdf) of distribution from mathematica
    ###### this should only work with analytic approximation, as otherwise q can only be accessed by a table, offering discrete values
    ###### then no need to generate all this!
    ###### however, what is a CDF of a distribution with pole at 0? (sure, should work as long as int(p(x)dx)=1
    if computation['p_theory']:
        results['p_range'] = ncid.variables['p_range'][:]      # distribution x-values
        results['p_exact'] = ncid.variables['p_exact'][:]      # distribution y-values
        #results['p_approx'] = ncid.variables['p_approx'][:]      # distribution y-values
        #results['alpha'] = ncid.variables['alpha'][:]
        #results['sigma_V'] = ncid.variables['sigma_V'][:]
        #results['I'] = ncid.variables['I'][:]
        #results['chi'] = ncid.variables['chi'][:]
    if computation['p_theory_hist']:
        results['p_hist'] = ncid.variables['p_hist'][:]        # discretized distribution (as histogram)

    #results['poisson'] = ncid.variables['poisson'][:]
    #if computation['draw_from_theory']:
        #results['p_est_inf'] = ncid.variables['p_est_inf'][:]
        #results['p_k'] = ncid.variables['p_k'][:]

    if computation['draw_finite_time']:
        #results['p_est_T'] = ncid.variables['p_est_T'][:]
        results['p_bayes_est'] = ncid.variables['p_bayes_est'][:]
        results['N_AP'] = ncid.variables['N_AP'][:]
        results['KS'] = ncid.variables['KS'][:]
        results['KL'] = ncid.variables['KL'][:]

    ncid.close()

    #results['steps'] = len(results['p_range'])
    #print "steps: ", results['steps']


    ### bootstrapping
    #results['bs_bayes'] = np.zeros((results['steps'],3))
    #for i in range(results['steps']):
        #results['bs_bayes'][i] = bootstrap(results['p_bayes_est'][:,:,i],100)


    #bayes = np.nansum(results['p_bayes_est'],axis=0)/(d_nu*N)
    ### plot the stuff
    if plot:
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        mpl.rcParams['font.size'] = 12
        mpl.rcParams['xtick.labelsize'] = 12
        mpl.rcParams['ytick.labelsize'] = 12

        if plot == 1:

            print("alpha: ", alpha_0)
            results['alpha'] = np.sqrt(results['alpha']**2+alpha_0**2)

            results['I_range'] = np.linspace(-0.5,0.2,801)
            results['f_I'] = 1./(2*math.pi*np.sqrt(model['tau_G']*model['tau_M']))*np.exp(-(results['I_range'])**2/(2*results['sigma_V']**2))
            results['I_distr'] = 1./np.sqrt(2*math.pi*results['alpha']**2)*np.exp(-(results['I_range']-results['I'])**2/(2*results['alpha']**2))
            results['tau_M_factor'] = 0.01/model['tau_M']
            return results

        if plot == 2:


            plt.figure(figsize=(10,6))
            ax1 = plt.axes([0.15,0.475,0.75,0.5])
            ax2 = ax1.twinx()
            ax3 = plt.axes([0.15,0.12,0.75,0.35])
            ax4 = plt.axes([0.65,0.25,0.21,0.17])


            max_1 = min(np.nanmax(results['p_exact']),5)*1.1

            #if computation['draw_finite_time']:
                #ax2.hist(results['rate_T'][0,0],bins=results['hist_range'],color='r',alpha=0.5,label=r'$\displaystyle \nu_{T=%g}^{N=%d}$'%(T,N))

            for k in range(computation['draw_from_theory']):

                #if computation['draw_from_theory']:
                    #ax1.plot(results['p_range'],results['p_est_inf'][k],'r',linewidth=2)

                for j in range(computation['draw_finite_time']):

                    if computation['draw_finite_time']:
                        #ax1.plot(results['p_range'],results['p_est_T'],'r--',linewidth=2)
                        ax1.plot(results['p_range'],results['p_bayes_est'][k,j],'k--',linewidth=3,label=r'$\displaystyle \rho(\nu)$ bayes')

                        #results['rate_inf_hist'] = np.histogram(results['rate_inf'][k],bins=results['hist_range'])
                        #ax2.hist(results['rate_inf'],bins=results['hist_range'],color='k',alpha=0.5,label=r'$\displaystyle \nu_{T\rightarrow\infty}^{N=%d}$'%N)
                        #ax2.plot(results['hist_range'][:-1],results['rate_inf_hist'][0],'k',drawstyle='steps-post',linewidth=2)
                        #ax2.plot(results['rate_inf'],np.ones(len(results['rate_inf'])),'ro')

                        #results['rate_inf_err'] = results['rate_inf_hist'][0]-results['p_hist']
                        #ax3.plot(results['hist_range'][:-1],results['rate_inf_err'],'k',drawstyle='steps-post')

                        #ax4.hist(results['rate_inf_err'],bins=np.linspace(-0.1*N,0.1*N,0.2*N+1),color='k',alpha=0.5)


                        #ax1.plot(results['p_range'],results['p_est_T'],'k--',linewidth=5)

                        #results['rate_T_hist'] = np.histogram(results['rate_T'],bins=results['hist_range'])


                        #ax2.plot(results['hist_range'][:-1],results['rate_T_hist'][0],'k',drawstyle='steps-post',linewidth=2)

                        #results['rate_T_err'] = results['rate_T_hist'][0]-results['p_hist']
                        #ax3.plot(results['hist_range'][:-1],results['rate_T_err'],'r',drawstyle='steps-post')

                        #ax4.hist(results['rate_T_err'],bins=np.linspace(-0.1*N,0.1*N,0.2*N+1),color='r',alpha=0.5)

            ax1.plot(results['p_range'],results['p_exact'],'k',linewidth=3,label=r'$\displaystyle \rho(\nu)$ theory')

            #ax2.step(results['hist_range'][:-1],results['p_hist'],color='k',where='post',linewidth=3,label=r'$\displaystyle \nu_{T=%g}^{N=%d}$'%(T,N))
            ax2.spines['top'].set_color('none')
            ax2.yaxis.set_ticks_position('right')
            ax2.yaxis.set_label_position('right')
            ax2.set_ylabel(r'bin counts H')
            ax2.set_xlim([0,5])
            ax2.set_ylim([0,l_H*N*max_1]) # this rescales it such that pdf and histogram are displaying the same scale
            ax1.legend(prop={'size':12})
            ax2.legend(prop={'size':12})

            ax1.set_ylim([0,max_1])
            #ax1.spines['right'].set_color('none')
            ax1.spines['top'].set_color('none')

            ax1.yaxis.set_ticks_position('left')
            ax1.yaxis.set_label_position('left')

            ax1.set_ylabel(r'$\displaystyle \rho(\nu)$')
            #ax1.set_xlabel(r'$\displaystyle \nu$ (Hz)')
            ax1.set_xticks([])

            ax1.set_xlim([0,5*rate])

            ax3.set_xlim([0,5*rate])
            ax3.set_ylim([-0.1*N,0.1*N])

            ax4.set_ylim([0,20])
            ax4.set_xlabel(r'$\displaystyle \Delta H$')

            #if save:
                #sv_name = './../pics/draw_samples_T=%g_N=%d_alpha=%g_eval.pdf' % (T,N,alpha_0)
                #plt.savefig(sv_name)
                #print "Figure saved as %s" % sv_name
            plt.show(block=False)

            plt.figure()
            plt.plot(results['p_range'][1:],np.cumsum(results['p_exact'][1:]*1/T),'k')
            #if computation['draw_from_theory']:
                #plt.plot(results['p_range'],np.cumsum(results['p_est_inf']*1/T),'r')
            if computation['draw_finite_time']:
                for k in range(computation['draw_from_theory']):
                    for j in range(computation['draw_finite_time']):
                        plt.plot(results['p_range'],np.cumsum(results['p_bayes_est'][k,j]*1/T),'r--')
            plt.show(block=False)

            if computation['draw_finite_time']:

                for k in range(computation['draw_from_theory']):
                    for j in range(computation['draw_finite_time']):
                        #print "K-S-test: ", np.max(np.abs(np.cumsum(results['p_exact'][1:]*1/T)-np.cumsum(results['p_bayes_est'][k,j,1:]*1/T)))
                        KL_log = np.log2(results['p_exact'][1:]/results['p_bayes_est'][k,j,1:])
                        mask = ~np.isinf(KL_log)
                        #print np.log((p_exact/T)/p_est_inf)
                        print("KL-Entropie: ", np.sum(results['p_exact'][1:][mask]*KL_log[mask]), "c++: ", results['KL'][k][j])


            #plt.figure()
            #plt.plot(results['p_range'],results['p_exact'],'k')
            #plt_bootstrap(results['p_range'],results['bs_bayes'],plt,'r')
            #plt.show(block=False)


                plt.figure()
                plt.plot(np.arange(computation['draw_from_theory']),results['KS'],'ko')
                plt.plot(np.arange(computation['draw_from_theory']),results['KL'],'ro')
                plt.xlim([-0.5,computation['draw_from_theory']-0.5])
                plt.ylim([0,0.5])
                plt.show(block=False)

            #plt.figure()
            #plt.xlim([-0.5,computation['draw_from_theory']-0.5])
            #plt.show(block=False)

            #poisson_range = np.linspace(0,5,T*5*5)
            #print poisson_range.shape
            #print results['poisson'].shape

            #plt.figure()
            #plt.plot(poisson_range,results['poisson'])
            #plt.show()

        #return results










def bootstrap(X,N,err=False,lower=0.025,upper=0.975,n=None):

  assert len(X.shape) <= 2, "Please only use 1/2-dim arrays (or 1 dim with error) for bootstrapping"

  if len(X.shape) == 2:
      X = np.array([x for sub in X for x in sub])


  if len(X.shape) == 1:
    mask = np.invert(np.isnan(X))
    if not np.sum(mask):
      return [np.nan]*3
    X_masked = X[mask]

    bootstrap_distr = np.zeros(N)

    for i in range(N):
      bootstrap_distr[i] = np.mean(bootstrap_resample(X_masked))

  elif len(X.shape) == 2:
    mask_nan = np.invert(np.isnan(X[0]))
    mask_inf = np.invert(np.isinf(X[1]))
    if not np.sum(mask_nan&mask_inf):
      return np.array([np.nan,np.nan,np.nan])

    X_masked = X[:,mask_nan&mask_inf]

    bootstrap_distr = np.zeros(N)

    for i in range(N):
      resample_X = bootstrap_resample(X_masked)
      bootstrap_distr[i] = np.average(resample_X[0],weights=1./resample_X[1])

  bootstrap_distr = np.sort(bootstrap_distr)

  #plt.figure()
  #plt.hist(bootstrap_distr,bins=100)
  #plt.show(block=False)

  lower_bound = bootstrap_distr[int(lower*N)]
  upper_bound = bootstrap_distr[int(upper*N)]
  mean = np.mean(bootstrap_distr)

  return np.array([mean, lower_bound, upper_bound])


def bootstrap_resample(X,n=None):

  if n == None:
    n = X.shape[-1]

  resample_i = np.floor(np.random.randint(n,size=n)).astype(int)

  if len(X.shape) == 1:
    X_resample = X[resample_i]

  elif len(X.shape) == 2:
    X_resample = X[:,resample_i]

  return X_resample


def plt_bootstrap(X,Y,ax,col,ls='-',ms='o',label=None,fit_func=None,p0=None,mask=None):

  assert len(X.shape) == 1, 'Please specify a one dimensional range array!'
  assert Y.shape[0] == X.shape[0], 'X and Y arrays should have the same length!'
  assert Y.shape[1] == 3, 'Y array should include mean value and upper and lower bound of confidence intervall!'

  if mask == None:
    mask = np.ones(len(X)).astype('bool')

  ax.plot(X[mask],Y[mask,0],'-',color=col,linestyle=ls,marker=ms,label=label,linewidth=2)

  ax.fill_between(X[mask],Y[mask,1],Y[mask,2],color=col,linestyle=ls,alpha=0.2,edgecolor=None,linewidth=0)

  if fit_func:
    Y_sigma = np.max(Y[:,1:] - Y[:,0].reshape(len(Y),1),axis=1)

    popt,pcov = curve_fit(fit_func,X[mask],Y[mask,0],sigma=Y_sigma[mask],p0=p0)

    perr = np.sqrt(np.diag(pcov))

    print('fit results: ', popt)
    print('fit errors: ', perr)
    ax.plot(X,fit_func(X,popt[0],popt[1]),'--',color='r',linewidth=0.5)

    return popt,perr







    #gamma = np.zeros((steps,steps,2))
    #chi = np.zeros((steps,steps,2))

    #gamma[:] = ncid.variables['gamma'][:].transpose(ax_list)
    #chi[:] = ncid.variables['chi'][:].transpose(ax_list)
    #gamma = gamma[:,:,0]
    #chi = chi[:,:,0]

    #plot_matrix = np.zeros((steps,steps))

    #mask_inconsistent = (gamma == -1)
    ##mask_no_peak = (gamma == -2)
    #mask_dark_matter = (gamma < 1)
    #plot_gamma = masked_array(1 - gamma**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
    #plot_chi = masked_array(chi,mask_inconsistent + mask_no_peak + mask_dark_matter)
    #plot_regions = masked_array(gamma,np.invert(mask_inconsistent + mask_no_peak))

    #ncid.close()
