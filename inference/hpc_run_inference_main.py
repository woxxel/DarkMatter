import sys, time, pickle
from pathlib import Path
from inference.BayesModel import *
from empirical.readData import *
import ultranest


def distribution_tests(BM,ref_values,results,plotting=False):
    ## compute KL-divergence and KS-test

    N_max = 100*BM.T

    N_array = np.arange(0,N_max)

    # res_values
    params_dict = {}
    for i,key in enumerate(BM.paramNames):
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




def run_inference(gamma, delta, nu_max, T, N, pathResults='.data/results/'):
    t_start = time.time()

    ## generate artificial data
    ref_values = {
        'gamma': gamma,
        'delta': delta,
        'nu_max': nu_max
    }
    nAnimals = 5

    mP = ModelParams(mode='artificial',parameter=ref_values,T=T,N=N,nAnimals=nAnimals)
    assert not (mP.rates is False), "No rates generated for given parameters"

    ## initialize Bayes model for inference
    BM = BayesModel(mP)
    BM.set_logLevel(logging.ERROR)
    BM.prepare_data(mP,mode='rates')
    BM.set_priors(hierarchical=[],two_pop=False)

    my_prior_transform = BM.set_prior_transform(vectorized=True)
    my_likelihood = BM.set_logl(vectorized=True,withZeros=True)

    sampler = ultranest.ReactiveNestedSampler(
        BM.paramNames, 
        my_likelihood, my_prior_transform,
        wrapped_params=BM.wrap,
        vectorized=True,num_bootstraps=20,
        ndraw_min=512
    )

    logger = logging.getLogger('ultranest')
    logger.setLevel(logging.ERROR)

    num_samples = 100

    sampling_result = sampler.run(
        min_num_live_points=num_samples,
        max_iters=20000,cluster_num_live_points=20,max_num_improvement_loops=3,
        show_status=True,viz_callback=False)

    ## calculating KL and KS divergence
    KL, KS = distribution_tests(BM,ref_values,sampling_result)

    res = {
        'empirical': {
            'parameters': ref_values,
            'mean_rates': mP.rates.mean(axis=0),
            'max_rates': mP.rates.max(axis=0),
        },
        'mean': sampling_result['posterior']['mean'],
        'stdev': sampling_result['posterior']['stdev'],
        'KL': KL,
        'KS': KS,
        'time_taken': time.time()-t_start
    }

    ## create file name, such that parameters are included and merely repeat numbers are added
    fileName = str(Path(pathResults) / 'test_inference_results')

    for var in [gamma,delta,nu_max,T,N]:
        fileName += f'_{var=:d}' if isinstance(var,int) else f'_{var=:.2f}'
    print(fileName)
    i=0
    while Path(f'{fileName}_{i}.pkl').is_file():
        i+=1

    pickle.dump(res,open(f'{fileName}_{i}.pkl','wb'))


if __name__ == '__main__':

    print('ID:',os.environ['SLURM_ARRAY_TASK_ID'])

    _, gamma, delta, nu_max, T, N = sys.argv
    run_inference(float(gamma), float(delta), float(nu_max), float(T), float(N))

