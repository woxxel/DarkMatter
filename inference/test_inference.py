from inference.utils.utils import adaptive_integration, p_nu, f
from inference.BayesModel import *
from empirical.readData import *

import pickle



def test_inference(n_repeat=5,ref_values=None,iter_key='N',iter_vals=[10,20,50,100,200,500,1000]):

    '''
        for each set of data, create artificial data, run inference on it and store relevant information (which one is that?)

        parameters to test are:
            - gamma
            - delta
            - nu_max
            - T
            - N
        how do I decide the remainder parameters?
    '''
    
    if ref_values is None:
        ref_values = {
            'gamma': 1.7,
            'delta': 4.5,
            'nu_max': 25.,
            'T': 1200,
            'N': 1000,
            'nAnimals': 5
        }
    ## test N

    # res = {}
    for val in iter_vals:
        
        ref_values[iter_key] = val

        # res[val] = {
        #     'mean': np.zeros((n_repeat,3)),
        #     'stdev': np.zeros((n_repeat,3)),
        #     'KL': np.zeros(n_repeat),
        #     'KS': np.zeros(n_repeat),
        #     'ref': ref_values
        # }


        for n in range(n_repeat):
            mP = ModelParams(mode='artificial',parameter=ref_values,T=ref_values['T'],N=ref_values['N'],nAnimals=ref_values['nAnimals'])
            if mP.rates is False:
                continue

            BM = BayesModel(mP)
            BM.set_logLevel(logging.ERROR)
            # BM.prepare_data(mP,mode='rates',key='WT')
            BM.prepare_data(mP,mode='rates')
            BM.set_priors(hierarchical=[],two_pop=False)

            my_prior_transform = BM.set_prior_transform()
            my_likelihood = BM.set_logl(withZeros=True)

            sampler = ultranest.ReactiveNestedSampler(
                BM.paraNames,
                my_likelihood, my_prior_transform,
                wrapped_params=BM.wrap,
                vectorized=True,num_bootstraps=20
            )

            logger = logging.getLogger("ultranest")
            logger.setLevel(logging.ERROR)

            num_samples = 400

            sampling_result = sampler.run(
                min_num_live_points=num_samples,
                max_iters=20000,cluster_num_live_points=20,max_num_improvement_loops=3,
                show_status=True,viz_callback=False)
            

            # res[val]['mean'][n,:] = sampling_result['posterior']['mean']
            # res[val]['stdev'][n,:] = sampling_result['posterior']['stdev']
            res = {}
            res['mean'] = sampling_result['posterior']['mean']
            res['stdev'] = sampling_result['posterior']['stdev']

            # res[val]['KL'][n], res[val]['KS'][n] = distribution_tests(BM,ref_values,sampling_result)
            res['KL'], res['KS'] = distribution_tests(BM,ref_values,sampling_result)

            pickle.dump(res,open(f'./results/test_inference_results_{iter_key}={val}_{n}.pkl','wb'))

            ## add additional parameter: duration for each run
            ## number of animals "pooled"
            ## maybe also: store maximum value of rates (and mean?)

            ## maybe adjust drawing samples: compute integral differently?!

            ## store cdf / KS-difference along with it?!

    return res

def distribution_tests(BM,ref_values,results,plotting=False):
    ## compute KL-divergence and KS-test

    N_max = 100*BM.T

    N_array = np.arange(0,N_max)

    # res_values
    params_dict = {}
    for i,key in enumerate(BM.paraNames):
        params_dict[key] = results['posterior']['mean'][i]

    p_nu_res = p_nu((N_array+0.5)/BM.T,params_dict)
    p_nu_res[~np.isfinite(p_nu_res)] = np.nan
    cum_res = np.nancumsum(p_nu_res)/np.nansum(p_nu_res)
    
    p_nu_ref = p_nu((N_array+0.5)/BM.T,ref_values)
    p_nu_ref[~np.isfinite(p_nu_ref)] = np.nan
    cum_ref = np.nancumsum(p_nu_ref)/np.nansum(p_nu_ref)
    
    if plotting:
        fig,ax = plt.subplots(2,1)
        ax[0].plot(N_array/BM.T,p_nu_ref,'k-',linewidth=2.5)
        ax[0].plot(N_array/BM.T,p_nu_res,'r--')

        ax[1].plot(N_array/BM.T,cum_ref,'k-',linewidth=2.5)
        ax[1].plot(N_array/BM.T,cum_res,'r--')

        plt.setp(ax[0],ylabel='probability')
        plt.setp(ax[1],xlim=ax[0].get_xlim(),xlabel='rate [Hz]',ylabel='cumulative probability')

        plt.show(block=False)

    KL = np.nansum(p_nu_res*np.log(p_nu_res/p_nu_ref))

    KS = np.nanmax(np.abs(cum_res-cum_ref))

    return KL, KS


def plot_results(iter_key,res,log=False):

    fig, ax = plt.subplots(3,1,figsize=(5,5))
    
    # N_array = np.array(list(res.keys()))
    ref_vals = {
        'gamma': 1.7,
        'delta': 4.5,
        'nu_max': 25.
    }

    n_repeat = None

    title_vals = {
        'gamma': r'$\gamma$',
        'delta': r'$\delta$',
        'nu_max': r'$\nu_{\max}$',
        'T': 'T',
        'N': 'N',
        'nAnimals': '# animals'
    }
    min_val = np.inf
    max_val = 0

    for j,var in enumerate(['gamma','delta','nu_max']):
        
        # ax[j].plot([0,0],[10,10],'k--')
        ax[j].axhline(ref_vals[var],color='k',linestyle='--')
        for i,val in enumerate(res.keys()):

            min_val = val if val<min_val else min_val
            max_val = val if val>max_val else max_val

            n_repeat = res[val]['mean'].shape[0]

            if log:
                ax[j].errorbar(val*(1+np.random.randn(n_repeat)/20.),res[val]['mean'][:,j],yerr=res[val]['stdev'][:,j],fmt='o')
            else:
                print('here')
                ax[j].errorbar(val+np.random.randn(n_repeat)/100.,res[val]['mean'][:,j],yerr=res[val]['stdev'][:,j],fmt='o')


        if log:
            ax[j].set_xscale('log')
        ax[j].spines[['top','right']].set_visible(False)
        plt.setp(ax[j],xlabel=title_vals[iter_key],ylabel=title_vals[var])
    plt.tight_layout()
    plt.show()
    


def plot_paper():

    '''
        define the figure as should be plotted for the paper
    '''
    steps = 20
    params = {
        'gamma': 1.7,
        'delta': 4.5,
        'nu_max': 25.,
        'T': 1200,
        'N': 100,
        'nAnimals': 5
    }
    
    fig = plt.figure(figsize=(12,8))

    ax_logp_points = fig.add_subplot(631)
    ax_logp_points_error = fig.add_subplot(634)
    ax_logp_max = fig.add_subplot(632)
    ax_logp_mean = fig.add_subplot(633)
    

    N_max = 100*params['T']
    nu_array = np.linspace(0,N_max/params['T'],N_max+1)


    ## count based
    N_max_steps = np.ceil(2*T/steps)*steps
    rateWnt_array = np.linspace(0,N_max_steps/params['T'],steps+1)
    p_N_AP = np.zeros(steps+1)
    for r,rateWnt in enumerate(rateWnt_array):

        p_N_AP[r] = adaptive_integration(f,0,params['nu_max'],
            args=(params,np.array([rateWnt*params['T']]),params['T'],True)
        )
    
    offset = 0.5
    p_N_AP_cont = p_nu(nu_array + offset,params,two_pop=True) / params['T']
    p_N_AP_cont_points = p_nu(rateWnt_array + offset,params,two_pop=True) / params['T']

    ax_logp_points.plot(nu_array,p_N_AP_cont,'k')
    ax_logp_points.plot(rateWnt_array,p_N_AP_cont_points,'rx')
