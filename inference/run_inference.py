import sys, time, pickle
from pathlib import Path
from inference.BayesModel import *
from empirical.readData import *

from scipy import interpolate, ndimage


def run_inference(ref_values,
        correct_N=5,bias_to_expected_max=0,bias_to_mean=0,
        path_results='data/results/',suffix='',
        sample_mode = 'dynesty',nP=8):

    t_start = time.time()

    print('in path (run_inference):',path_results)

    T = ref_values['T']
    N = ref_values['N']
    nAnimals = 5

    mP = ModelParams(mode='artificial',path_results=path_results,parameter=ref_values,T=T,N=N,nAnimals=nAnimals,suffix=suffix)
    assert not (mP.rates is False), "No rates generated for given parameters"

    ## initialize Bayes model for inference
    BM = BayesModel(mP)
    BM.set_logLevel(logging.ERROR)
    BM.prepare_data(mP,mode='rates')
    BM.set_priors(hierarchical=[],two_pop=mP.two_pop)


    my_prior_transform = BM.set_prior_transform(vectorized=sample_mode=='ultranest')
    my_likelihood = BM.set_logl(vectorized=sample_mode=='ultranest',
        correct_N=correct_N,
        bias_to_expected_max=bias_to_expected_max,
        bias_to_mean=bias_to_mean,
    )

    result_values = {
        'mean': {},
        'CI': {}
    }

    num_samples = 100
    if sample_mode=='ultranest':

        import ultranest
        from ultranest.popstepsampler import PopulationSliceSampler, generate_region_oriented_direction
        
        sampler = ultranest.ReactiveNestedSampler(
            BM.paramNames, 
            my_likelihood, my_prior_transform,
            wrapped_params=BM.wrap,
            vectorized=True,num_bootstraps=20,
            ndraw_min=512
        )

        logger = logging.getLogger('ultranest')
        logger.setLevel(logging.ERROR)

        sampler.stepsampler = PopulationSliceSampler(popsize=20,nsteps=10,
                generate_direction=generate_region_oriented_direction)

        sampling_result = sampler.run(
            min_num_live_points=num_samples,
            max_iters=20000,cluster_num_live_points=20,max_num_improvement_loops=1,
            show_status=True,viz_callback=False
        )

        ## calculating KL and KS divergence
        
        for i,key in enumerate(BM.paramNames):
            result_values['mean'][key] = sampling_result['posterior']['mean'][i]
            result_values['CI'][key] = [
                sampling_result['posterior']['mean'][i] - sampling_result['posterior']['stdev'][i],
                sampling_result['posterior']['mean'][i] + sampling_result['posterior']['stdev'][i],
            ]
        
    elif sample_mode=='dynesty':

        from dynesty import NestedSampler, pool as dypool, utils as dyfunc

        if nP>1:
            with dypool.Pool(nP,my_likelihood,my_prior_transform) as pool:
                
                sampler = NestedSampler(pool.loglike,pool.prior_transform,BM.nParams,
                        pool=pool,
                        nlive=num_samples,
                        bound='single',
                        reflective=[BM.priors['p']['idx']] if mP.two_pop else False,
                        sample='rslice'
                    )
                sampler.run_nested(dlogz=1.)
        else:

            sampler = NestedSampler(my_likelihood,my_prior_transform,BM.nParams,
                nlive=num_samples,
                bound='single',
                reflective=[BM.priors['p']['idx']] if mP.two_pop else False,
                sample='rslice'
            )
            sampler.run_nested(dlogz=1.)
        sampling_result = sampler.results

        # samples,weights = sampler.results.samples, sampler.results.importance_weights()

        posterior = build_posterior(BM,sampler.results)


        # mean,cov = dyfunc.mean_and_cov(samples,weights)

        for i,key in enumerate(BM.paramNames):
            result_values['mean'][key] = posterior[key]['mean']
            result_values['CI'][key] = posterior[key]['CI'][[1,3]]

    if mP.two_pop:
        for key in ['mean','CI']:
            result_values[key]['nu_max_2'] = result_values[key]['nu_max_1']

    KL, KS = distribution_tests(BM,ref_values,result_values['mean'])

    '''
        further, store 
            * number of silent neurons
            * posterior distribution of parameters
    '''

    res = {
        'empirical': {
            'parameters': ref_values,
            'mean_rates': mP.rates.mean(axis=0),
            'max_rates': mP.rates.max(axis=0),
            'silent_neurons': (mP.rates==0).sum(axis=0),
            'rates': mP.rates,
        },
        'posterior': posterior,
        'KL': KL,
        'KS': KS,
        'time_taken': time.time()-t_start
    }

    ## create file name, such that parameters are included and merely repeat numbers are added
    if path_results:
        fileName = str(Path(path_results,'data','inference_test'))

        mode_str = '_modes'
        if correct_N>0:
            mode_str += f'_{correct_N=}'
        if bias_to_expected_max>0:
            mode_str += f'_{bias_to_expected_max=}'
        if bias_to_mean>0:
            mode_str += f'_{bias_to_mean=}'
        
        fileName += mode_str

        for key,val in ref_values.items():
            fileName += f'_{key}={val:d}' if isinstance(val,int) else f'_{key}={val:g}'
        i=0
        while Path(f'{fileName}_{i}.pkl').is_file():
            i+=1
        print('saving to:',f'{fileName}_{i}.pkl')

        pickle.dump(res,open(f'{fileName}_{i}.pkl','wb'))
    else:
        return BM, sampling_result, sampler, res


def build_posterior(BM,results,nsteps=101,smooth_sigma=1,plot=False):

    posterior = {}
    
    for i,key in enumerate(BM.paramNames):
        samp = results.samples[:,i]

        # sort samples
        samples_sorted = np.sort(samp)
        idx_sorted = np.argsort(samp)


        # get corresponding weights
        weights = results.importance_weights()
        sw = weights[idx_sorted]

        cumsw = np.cumsum(sw)

        quants = np.interp([0.001,0.05,0.341,0.5,0.841,0.95,0.999], cumsw, samples_sorted)
        
        low,high = quants[[0,-1]]
        x = np.linspace(low,high,nsteps)

        f = interpolate.interp1d(samples_sorted,cumsw,bounds_error=False,fill_value='extrapolate')
        
        posterior[key] = {
            'CI': quants[1:-1],
            'mean': quants[3],
            'x': x[:-1],
            'p_x': ndimage.gaussian_filter1d(f(x[1:]) - f(x[:-1]),smooth_sigma) if smooth_sigma>0 else f(x[1:]) - f(x[:-1])
        }

    return posterior


def plot_posterior(posterior):

    two_pop = len(posterior.keys())>3
    fig,ax = plt.subplots(3,2 if two_pop else 1,figsize=(12,8))
    
    for i,key in enumerate(posterior.keys()):
        axx = ax[i//2][i%2] if two_pop else ax[i]
        for j,CI_val in enumerate(posterior[key]['CI']):
            axx.axvline(CI_val,color='g' if j==2 else 'r',linestyle='-' if j==2 else '--')
        axx.plot(posterior[key]['x'],posterior[key]['p_x'])

        
        axx.set_title(f"{key} = {posterior[key]['mean']:.2f} ({posterior[key]['CI'][1]:.2f} , {posterior[key]['CI'][3]:.2f})")
        axx.spines[['top','right']].set_visible(False)
        plt.setp(axx,yticks=[],xlabel=key)

    plt.tight_layout()
    plt.show()



def distribution_tests(BM,ref_values,result_values,plotting=False):
    ## compute KL-divergence and KS-test

    N_max = 100*BM.T

    N_array = np.arange(0,N_max)

    # res_values
    offset = 0.5
    nu_array = (N_array+offset)/BM.T

    p_nu_res = get_p_nu(nu_array,result_values,BM.T,correct_N=10)
    cum_res = np.nancumsum(p_nu_res)/np.nansum(p_nu_res)
    
    p_nu_ref = get_p_nu(nu_array,ref_values,BM.T,correct_N=10)
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

    KS_all = np.nanmax(np.abs(cum_res-cum_ref))
    
    P = lambda nu : p_nu(nu,ref_values)
    Q = lambda nu : p_nu(nu,result_values)

    KL_all = adaptive_integration(lambda nu: P(nu)*np.log(P(nu)/Q(nu)),0,min(ref_values['nu_max_1'],result_values['nu_max_1']),args=(),eps_pow=-8,eps_thr=-4)

    ## get empirical ecdf
    plotting = False
    if plotting:
        fig,ax = plt.subplots(5,1,figsize=(5,12))
        for a in range(BM.rates.shape[1]):
            x = np.sort(BM.rates[:,a])
            y = np.arange(len(x))/float(len(x))

            ax[a].plot(x,y,'r-')

            pdf_theory = np.zeros(BM.rates.shape[0])
            ecdf_theory = np.zeros(BM.rates.shape[0])
            for i,nu in enumerate(x):
                idx = np.argmin(np.abs(nu-nu_array))
                # print(i,idx)
                pdf_theory[i] = p_nu_res[idx]
                ecdf_theory[i] = cum_res[idx]
            ax[a].plot(x,ecdf_theory,'k--')

            KS = np.nanmax(np.abs(y-ecdf_theory))
            print(f'{a=}: {KS=} vs {KS_all=}')

            # print(pdf_theory)
            KL = np.nansum(pdf_theory[1:]*np.log(pdf_theory[1:]/y[1:]))
            print(f'{a=}: {KL=} vs {KL_all=}')

            # KL = adaptive_integration(lambda nu: P(nu)*np.log(P(nu)/Q(nu)),0,min(ref_values['nu_max'],result_values['nu_max']),args=(),eps_pow=-8,eps_thr=-4)

        plt.show(block=False)


    # KL = np.nansum(p_nu_res*np.log(p_nu_res/p_nu_ref))
    # print(f'{KL_1=}, {KL=}')

    return KL_all, KS_all


if __name__ == '__main__':

    # print('ID:',os.environ['SLURM_ARRAY_TASK_ID'])
    suffix = f"{os.environ['SLURM_ARRAY_JOB_ID']}_{os.environ['SLURM_ARRAY_TASK_ID']}"
    print('suffix:',suffix)

    parameter_string, path_results, correct_N,bias_to_expected_max, bias_to_mean = sys.argv[1:]
        
    pairs = parameter_string.split('__')
    ref_values = {}
    for pair in pairs:
        key,val = pair.split('=')
        ref_values[key] = float(val)
    
    print('input parameters: ',ref_values)

    run_inference(ref_values,
        correct_N=int(correct_N),bias_to_expected_max=float(bias_to_expected_max),bias_to_mean=float(bias_to_mean),
        path_results=path_results,suffix=suffix)