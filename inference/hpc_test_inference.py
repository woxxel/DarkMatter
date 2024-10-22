from pathlib import Path
# import time, inspect
from inference.utils.utils import adaptive_integration, p_nu, f

import pickle
from .utils.connections import *



def test_inference(n_repeat=5,ref_values=None,iter_key='N',iter_vals=[10,20,50,100,200,500,1000],ssh_config_file_name='id_ed25519_GWDG'):

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
    
    nAnimals = 5
    if ref_values is None:
        ref_values = {
            'gamma': 1.7,
            'delta': 4.5,
            'nu_max': 25.,
            'T': 1200,
            'N': 1000,
        }
    
    cpus = 8
    # mouse = os.path.basename(path_source)

    # if suffix:
    #     suffix = f'_{suffix}' if not suffix.startswith('_') else suffix
    path_results = '/scratch/users/schmidt124/DarkMatter/inference'
    path_code = '~/program_code/DarkMatter'
    submit_file = f"{path_code}/sbatch_submit.sh"

    # if path_target is None:
    #     path_target = path_source

    ## setting up connection to server
    username = 'schmidt124'
    proxyServerName = 'login.gwdg.de'
    serverName = 'login-dbn02.hpc.gwdg.de'
    
    ssh_key_file = f'/home/wollex/.ssh/{ssh_config_file_name}'

    client = establish_connection(serverName,username,ssh_key_file,proxyJump=proxyServerName)
    sftp_client = client.open_sftp()

    inference_script = f'{path_code}/run_inference.py'

    ## first, write the neuron detection script to the server
    with open('inference/hpc_run_inference_main.py','r') as f_open:
        lines = f_open.readlines()
    _, stdout, stderr = client.exec_command(f"""cat > {inference_script} <<- EOF
{("").join(lines)}
    EOF
    """)
    
    ## iterate through specified parameters
    for val in iter_vals:
        ref_values[iter_key] = val

        # for _ in range(n_repeat):
            
#SBATCH -o {path_code}/logs/inference_%A_%a.log
#SBATCH -e {path_code}/logs/inference_%A_%a.log
        _, stdout, stderr = client.exec_command(f"""cat > {submit_file} <<- EOF
#!/bin/bash -l
#SBATCH -J distribution_inference
#SBATCH -a 1-{n_repeat}
#SBATCH -A cidbn_legacy
#SBATCH -p cidbn
#SBATCH -c {cpus}
#SBATCH -t 01:00:00
#SBATCH -o {path_results}/logs/inference_%A_%a.log
#SBATCH -e {path_results}/logs/inference_%A_%a.log
#SBATCH --mem=16000

# module load gsl
# module load netcdf-c

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS=1

python3 {inference_script} {ref_values['gamma']} {ref_values['delta']} {ref_values['nu_max']} {ref_values['T']} {ref_values['N']} {path_results}/data
EOF
""")
        _, stdout, stderr = client.exec_command(f'/usr/local/slurm/current/install/bin/sbatch {submit_file}')
        print(stdout.read(),stderr.read())


        ## maybe adjust drawing samples: compute integral differently?!

    # return res


def load_results(iter_key,iter_val,n_repeat):

    fileName = f'.data/results/test_inference_results_{iter_key}={iter_val}_{n_repeat}.pkl'
    if Path(fileName).is_file():
        with open(f'.data/results/test_inference_results_{iter_key}={iter_val}_{n_repeat}.pkl','rb') as f_open:
            return pickle.load(f_open)
    else:
        print('file not found')
    


def plot_results(iter_key,iter_vals,n_repeat,log=False):

    fig, ax = plt.subplots(3,1,figsize=(5,5))
    
    # N_array = np.array(list(res.keys()))
    ref_vals = {
        'gamma': 1.7,
        'delta': 4.5,
        'nu_max': 25.
    }

    title_vals = {
        'gamma': r'$\gamma$',
        'delta': r'$\delta$',
        'nu_max': r'$\nu_{\max}$',
        'T': 'T',
        'N': 'N',
        'nAnimals': '# animals'
    }
    # min_val = np.inf
    # max_val = 0

    for j,var in enumerate(['gamma','delta','nu_max']):
        
        # ax[j].plot([0,0],[10,10],'k--')
        ax[j].axhline(ref_vals[var],color='k',linestyle='--')

        for val in iter_vals:

            for n in range(n_repeat):
                res = load_results(iter_key,val,n)
                ax[j].errorbar(val+np.random.randn()/100.,res['mean'][j],yerr=res['stdev'][j],fmt='o')


        # for i,val in enumerate(res.keys()):

        #     # min_val = val if val<min_val else min_val
        #     # max_val = val if val>max_val else max_val

        #     n_repeat = res[val]['mean'].shape[0]

        #     if log:
        #         ax[j].errorbar(val*(1+np.random.randn(n_repeat)/20.),res[val]['mean'][:,j],yerr=res[val]['stdev'][:,j],fmt='o')
        #     else:
        #         print('here')
        #         ax[j].errorbar(val+np.random.randn(n_repeat)/100.,res[val]['mean'][:,j],yerr=res[val]['stdev'][:,j],fmt='o')


        if log:
            ax[j].set_xscale('log')
        ax[j].spines[['top','right']].set_visible(False)
        plt.setp(ax[j],xlabel=title_vals[iter_key],ylabel=title_vals[var])
    plt.setp(ax[0],ylim=[0.5,3])
    plt.setp(ax[1],ylim=[2,7])
    plt.setp(ax[2],ylim=[10,50])
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
