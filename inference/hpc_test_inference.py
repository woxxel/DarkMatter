from pathlib import Path
# import time, inspect
from inference.utils.utils import adaptive_integration, p_nu, f
import copy
import numpy as np
from matplotlib import pyplot as plt

import pickle
from .utils.connections import *

# import seaborn as sns

'''
    store posterior distribution along with results
'''

def test_inference(n_repeat=5,ref_values=None,
        iter_key='N',iter_vals=[10,20,50,100,200,500,1000],
        correct_N=5, bias_to_expected_max=0, bias_to_mean=0,
        ssh_config_file_name='id_ed25519_GWDG',
        **kwargs):

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
            'gamma_1': 1.7,
            'delta_1': 4.5,
            'nu_max_1': 25.,
            'T': 1200,
            'N': 1000,
        }

    for key in ref_values.keys():
        if key in kwargs.keys():
            ref_values[key] = kwargs[key]
    
    # print(ref_values)
    cpus = 8
    # mouse = os.path.basename(path_source)

    path_results = '/dbnscratch/users/schmidt124/DarkMatter/inference'
    path_code = '~/program_code/DarkMatter'
    submit_file = f"{path_code}/sbatch_submit.sh"


    ## setting up connection to server
    username = 'schmidt124'
    proxyServerName = 'login.gwdg.de'
    serverName = 'login-dbn02.hpc.gwdg.de'
    
    ssh_key_file = f'/home/wollex/.ssh/{ssh_config_file_name}'

    client = establish_connection(serverName,username,ssh_key_file,proxyJump=proxyServerName)

    ## hand over the inference script to the server
    inference_script = f'{path_code}/run_inference.py'
    with open(Path(__file__).parent / 'run_inference.py','r') as f_open:
        lines = f_open.readlines()
    _, stdout, stderr = client.exec_command(f"""cat > {inference_script} <<- EOF
{("").join(lines)}
EOF
""")
    
    ## iterate through specified parameters
    for val in iter_vals:
        ref_values[iter_key] = val

        print(ref_values)

        pass_parameter_string = '__'.join([f'{key}={val}' for key,val in ref_values.items()])

        # print(pass_parameter_string)
        # return
        _, stdout, stderr = client.exec_command(f"""cat > {submit_file} <<- EOF
#!/bin/bash -l
#SBATCH -J distribution_inference
#SBATCH -a 1-{n_repeat}%20
#SBATCH -A cidbn_legacy
#SBATCH -p cidbn
#SBATCH -c {cpus}
#SBATCH -t 10:00:00
#SBATCH -o {path_results}/logs/inference_%A_%a.log
#SBATCH -e {path_results}/logs/inference_%A_%a.log
#SBATCH --mem=16000

module use /usr/users/cidbn_sw/sw/modules
module load cidbn_dynesty-2.1.4_py-3.11

export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS=1

python3 {inference_script} {pass_parameter_string} {path_results} {correct_N} {bias_to_expected_max} {bias_to_mean}
EOF
""")
        _, stdout, stderr = client.exec_command(f'/usr/local/slurm/current/install/bin/sbatch {submit_file}')
        print(stdout.read(),stderr.read())



def load_results(base_path,parameters,n_repeat):

    fileName = base_path
    for key,val in parameters.items():
        fileName += f'_{key}={val:d}' if isinstance(val,int) else f'_{key}={val:g}'
    fileName = Path(f'{fileName}_{n_repeat}.pkl')

    # print(fileName)
    if Path(fileName).is_file():
        # print(fileName)
        # try:
        with open(fileName,'rb') as f_open:
            return pickle.load(f_open)
        # except:
            # print('error while loading')
    


def plot_results(iter_key,iter_vals,n_repeat,ref_values=None,
        correct_N=5, bias_to_expected_max=0, bias_to_mean=0,
        base_path='results_inference/inference_test_',save=False,log=False,**kwargs):
    

    mode_str = '_modes'
    if correct_N>0:
        mode_str += f'_{correct_N=}'
    if bias_to_expected_max>0:
        mode_str += f'_{bias_to_expected_max=}'
    if bias_to_mean>0:
        mode_str += f'_{bias_to_mean=}'

    base_path += mode_str
    # print(base_path)
    # N_array = np.array(list(res.keys()))
    if ref_values is None:
        ref_vals = {
            'gamma_1': 1.7,
            'delta_1': 4.5,
            'nu_max_1': 25.,
            'T': 1200,
            'N': 1000,
        }
    else:
        ref_vals = copy.deepcopy(ref_values)

    for key in ref_vals.keys():
        if key in kwargs.keys():
            ref_vals[key] = kwargs[key]
    print(ref_vals)

    para_keys = []
    for key in ref_vals.keys():
        if key!='N' and key!='T':
            para_keys.append(key)
    
    title_vals = {
        'gamma_1': r'$\gamma$',
        'delta_1': r'$\delta$',
        'nu_max_1': r'$\nu_{\max}$',
        'gamma_2': r'$\gamma_2$',
        'delta_2': r'$\delta$_2',
        'nu_max_2': r'$\nu_{\max}$',
        'p': r'$p$',
        'T': 'T',
        'N': 'N',
        'nAnimals': '# animals'
    }

    # print(para_keys)
    nParams = 4 if 'p' in para_keys else 3
    fig, ax = plt.subplots(4+nParams,1,figsize=(5,12),sharex=True)
    # print(nParams)
    ax_param_idx = 2
    ax_distr = ax[0]
    ax_finished = ax[1]

    ax_time = ax_finished.twinx()

    ax_KL = ax[ax_param_idx+nParams]
    ax_KS = ax[ax_param_idx+nParams+1]

    colors = ['tab:blue']
    if 'p' in para_keys:
        colors.append('tab:orange')

    x_range = (iter_vals[-1] - iter_vals[0])
    if log:
        x_offset = lambda x, factor : np.exp(factor * np.log(x))
        xlims = [iter_vals[0] / x_offset(x_range,0.1),iter_vals[-1] * x_offset(x_range,0.05)]
    else:
        x_offset = lambda x, factor : factor * x
        xlims = [iter_vals[0] - x_offset(x_range,0.25),iter_vals[-1] + x_offset(x_range,0.1)]
    prior_baseline = iter_vals[0] - x_offset(x_range,0.1)

    # one_distr = not iter_key.startswith(('gamma','delta','nu_max'))
    one_distr = False
    multi_distr = iter_key.startswith(('gamma','delta','nu_max'))
    if one_distr:

        nu_array = np.linspace(0,ref_vals['nu_max_1'],10001)
        rho = p_nu(nu_array,ref_vals)
        # if log:
            # rho = np.exp(rho)
        # ax_distr.plot(val - rho,nu_array,'k-',linewidth=0.5)
        ax_distr.fill_between(x=nu_array,y1=rho,y2=0,alpha=0.5,color='b')


    s_kwargs = {"s": 40, "marker": "_"}
    for i,val in enumerate(iter_vals):
    
        ref_vals[iter_key] = val

        if log:
            # val_offset =  
            # val_pos = val * np.exp(0.005*np.log(x_range))
            val_width = val * x_offset(x_range,0.01)
            val_scatter = val * np.exp( + np.random.randn(n_repeat)*0.005*np.log(x_range))
        else:
            val_width = x_offset(x_range,0.1)
            val_scatter = val + np.random.randn(n_repeat)*0.01*x_range
        
        v_kwargs = {"showextrema": False, "showmedians": False, "showmeans": True, "widths": val_width}
        
        ''' PLOT REFERENCE DISTRIBUTION'''
        if (multi_distr or i==0):
            nu_array = np.linspace(0,ref_vals['nu_max_1'],10001)
            rho = p_nu(nu_array,ref_vals)*x_offset(x_range,0.1)
            if log:
                rho = np.exp(rho)
            # ax_distr.plot(val - rho,nu_array,'k-',linewidth=0.5)
            ax_distr.fill_betweenx(y=nu_array,x1=val+rho,x2=val,alpha=0.5,color='b')

        # for idx,c in enumerate(colors):
        # for j,var in enumerate([f'gamma_{idx+1}',f'delta_{idx+1}',f'nu_max_{idx+1}','p']):
        for j,var in enumerate(para_keys):
            label_str = 'truth' if i==0 and j==0 else None
            if var==iter_key:
                ax[ax_param_idx+j%nParams].plot(iter_vals,iter_vals,'k--',linewidth=0.5,label=label_str)
            else:
                ax[ax_param_idx+j%nParams].axhline(ref_vals[var],color='k',linestyle='--',linewidth=0.5,label=label_str)
        

        result_blueprint = np.full(n_repeat,np.NaN)
        result_CI_blueprint = np.full((n_repeat,2),np.NaN)
        results = {
            'mean': {},
            # 'stdev': {},
            'CI': {},

            'KL': np.full(n_repeat,np.NaN),
            'KS': np.full(n_repeat,np.NaN),

            'time_taken': np.full(n_repeat,np.NaN)
        }
        for key in para_keys:
            results['mean'][key] = copy.deepcopy(result_blueprint)
            # results['stdev'][key] = copy.deepcopy(result_blueprint)
            results['CI'][key] = copy.deepcopy(result_CI_blueprint)

        # print(results['mean'].keys())
        for n in range(n_repeat):
            

            res = load_results(base_path,ref_vals,n)

            if res: 
                try:
                    for i,key in enumerate(para_keys):

                        if key.startswith('nu_max') and key[-1]!='1':#idx>0:
                            continue
                        # print(res['posterior'][key]['mean'])
                        results['mean'][key][n] = res['posterior'][key]['mean']
                        # results['stdev'][key][n] = res['posterior'][key]['CI']
                        results['CI'][key][n,:] = res['posterior'][key]['CI'][[1,-2]]
                except:
                    pass
                    # for key in para_keys:

                    #     results['mean'][key][n] = res['mean'][key]
                    #     results['stdev'][key][n] = res['stdev'][key]

                results['KL'][n] = res['KL']
                results['KS'][n] = res['KS']

                results['time_taken'][n] = res['time_taken']

        for idx,c in enumerate(colors):
            for j,var in enumerate([f'gamma_{idx+1}',f'delta_{idx+1}',f'nu_max_{idx+1}']):
            
                if np.all(np.isnan(results['mean'][var])):
                    continue
                
                if var.startswith('nu_max') and var[-1]!='1':#idx>0:
                    continue
                
                ax[ax_param_idx+j].errorbar(val_scatter,results['mean'][var],yerr=np.abs(results['mean'][var][:,np.newaxis] - results['CI'][var]).T,fmt='o',markersize=0.3,linewidth=0.2,color='gray')#,label='inference results' if (i==0 and j==0 and idx==0) else None)

                custom_violin(ax[ax_param_idx+j], results['mean'][var][np.isfinite(results['mean'][var])], val, c, c, 0.6, side="right", scatter_kwargs=s_kwargs, violin_kwargs=v_kwargs)

                if idx==0:
                    plt.setp(ax[ax_param_idx+j],ylabel=title_vals[var])
        
        if 'p' in para_keys and not np.all(np.isnan(results['mean']['p'])):
            var = 'p'
            c = 'tab:blue'
            # ax[ax_param_idx+3].errorbar(val_scatter,results['mean'][var],yerr=results['stdev'][var],fmt='o',markersize=0.3,linewidth=0.2,color='gray')
            ax[ax_param_idx+3].errorbar(val_scatter,results['mean'][var],yerr=np.abs(results['mean'][var][:,np.newaxis]-results['CI'][var]).T,fmt='o',markersize=0.3,linewidth=0.2,color='gray')
                    
            custom_violin(ax[ax_param_idx+3], results['mean'][var][np.isfinite(results['mean'][var])], val, c, c, 0.6, side="right", scatter_kwargs=s_kwargs, violin_kwargs=v_kwargs)

            plt.setp(ax[ax_param_idx+3],ylabel=title_vals[var])


        ''' PLOT TEST STATISTICS '''
        ax_KL.scatter(val_scatter,results['KL'],color='k',s=2)
        ax_KS.scatter(val_scatter,results['KS'],color='k',s=2)

        n_finished = np.isfinite(results['mean']['gamma_1']).sum()
        ax_finished.plot(val,n_finished,'ko')

        ax_time.scatter(val_scatter,results['time_taken'],color='grey',s=0.5)
    
    ax_finished.axhline(n_repeat,color='k',linestyle='--',linewidth=0.5)


    ''' PLOT PRIOR DISTRIBUTIONS '''
    exponential_fun = lambda x,mu,sigma : np.exp(-(x-mu)**2/(2*sigma**2))#/np.sqrt(2*np.pi*sigma**2)
    
    nu_array = np.linspace(0,100,1001)
    delta_array = np.linspace(0,10,1001)
    gamma_array = np.linspace(0.5,3,1001)

    gamma_prior = x_offset(x_range,0.1)*exponential_fun(gamma_array,2.,0.5)
    delta_prior = x_offset(x_range,0.1)*exponential_fun(delta_array,6.,2.)
    nu_prior = x_offset(x_range,0.1)*exponential_fun(nu_array,0.,50.)

    ax[ax_param_idx+0].fill_betweenx(y=gamma_array,x1=prior_baseline - gamma_prior,x2=prior_baseline,alpha=0.5,color='tab:green',label='prior')
    ax[ax_param_idx+1].fill_betweenx(y=delta_array,x1=prior_baseline - delta_prior,x2=prior_baseline,alpha=0.5,color='tab:green')
    ax[ax_param_idx+2].fill_betweenx(y=nu_array,x1=prior_baseline - nu_prior,x2=prior_baseline,alpha=0.5,color='tab:green')


    ''' ADJUST AXES PROPERTIES '''
    ax[ax_param_idx].legend(bbox_to_anchor=(.05, .95), loc='upper left', borderaxespad=0.)

    plt.setp(ax[ax_param_idx],ylim=[1.,3.5])
    plt.setp(ax[ax_param_idx+1],ylim=[2.5,10])
    plt.setp(ax[ax_param_idx+2],ylim=[10,60])
    plt.setp(ax_finished,ylim=[0,n_repeat*1.1],ylabel='finished runs')

    plt.setp(ax_KL,ylabel='KL div.')
    plt.setp(ax_KS,ylabel='KS stat.')

    if one_distr:
        plt.setp(ax_distr,xlim=[0,2],xlabel='rate [Hz]')
    else:
        plt.setp(ax_distr,ylim=[0,2],ylabel='rate [Hz]',xlim=xlims)

    ax_time.set_yscale('log')
    plt.setp(ax_time,ylim=[10**1,10**5],ylabel='computation time [s]',xlim=xlims)#,xticklabels=[])


    if log:# and not one_distr:
        ax_distr.set_xscale('log')

    for j,axx in enumerate(ax):
        axx.spines[['top','right']].set_visible(False)
        # plt.setp(axx,xlim=xlims)

        # if j<6:
            # plt.setp(axx,xticklabels=[])
        if j==(nParams+4-1):
            plt.setp(axx,xlabel=title_vals[iter_key])
    
    ax_time.spines[['top','left']].set_visible(False)
    
    plt.tight_layout(h_pad=0.1)

    if save:
        file_format='png'
        sv_name = f"./figures/inference_tests_pop={'2' if 'p' in para_keys else '1'}_{iter_key}.{file_format}"
        plt.savefig(sv_name)
        print('Figure saved as "%s"' % sv_name)
    # else:
    
    # plt.show(block=False)
    plt.show()
    


def plot_paper():

    '''
        define the figure as should be plotted for the paper
    '''
    steps = 20
    params = {
        'gamma_1': 1.7,
        'delta_1': 4.5,
        'nu_max_1': 25.,
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


def custom_violin(ax, data, pos, fc='b', ec='k', alpha=0.7, percentiles=[25, 50, 75], side="both", scatter_kwargs={}, violin_kwargs={}):
    """Customized violin plot.
    ax: axes.Axes, The axes to plot to
    data: 1D array like, The data to plot
    pos: float, The position on the x-axis where the violin should be plotted
    fc: color, The facecolor of the violin
    ec: color, The edgecolor of the violin
    alpha: float, The transparancy of the violin
    percentiles: array like, The percentiles to be marked on the violin
    side: string, Which side(s) of the violin should be cut off. Options: 'left', 'right', 'both'
    scatter_kwargs: dict, Keyword arguments for the scatterplot
    violin_kwargs: dict, Keyword arguments for the violinplot"""

    parts = ax.violinplot(data, positions=[pos], **violin_kwargs)
    for pc in parts['bodies']:
        m = np.mean(pc.get_paths()[0].vertices[:, 0])
        if side == "left":
            points_x = pos - 0.05
            pc.get_paths()[0].vertices[:, 0] = np.clip(pc.get_paths()[0].vertices[:, 0], -np.inf, m)
        elif side == "right":
            points_x = pos + 0.05
            pc.get_paths()[0].vertices[:, 0] = np.clip(pc.get_paths()[0].vertices[:, 0], m, np.inf)
        else:
            points_x = pos
        pc.set_facecolor(fc)
        pc.set_edgecolor(ec)
        pc.set_alpha(alpha)
        parts['cmeans'].set_color(fc)

    # perc = np.percentile(data, percentiles)
    # for p in perc:
    #     ax.scatter(points_x, p, color=ec, zorder=3, **scatter_kwargs)





def get_peak_position_and_height(params):

    # calculate the position of the distribution peak
    nu_peak = lambda gamma,delta,nu_max :  np.exp(np.log(nu_max) - ((gamma*delta)**2-2*((gamma)**2- 1) + gamma*delta*np.sqrt((gamma*delta)**2-4*((gamma)**2-1))) /(4*(gamma**2 - 1)**2))

    # get the maximum value of the distribution
    p_nu(nu_peak(**params),params)




#     if compile_it:
#         print('compiling...')
#         compile_path = Path('/home/uni02/UZDN/schmidt124/program_code/DarkMatter/theory')
#         compile_script = f"{path_code}/compile_cpp.sh"
#         _, stdout, stderr = client.exec_command(f"""
# cat > {compile_script} <<- EOF
# #!/bin/bash -l

# module load gsl
# module load netcdf-c
# echo "compiling generate_measures"

# g++ -g -Wall -ansi -o {compile_path / 'generate_measures'} {compile_path / 'generate_measures.cpp'} -I{compile_path / 'src' / '3rdParty' /'include'} -L{compile_path / 'src' / '3rdParty' / 'lib'}  -lgsl -lgslcblas -lnetcdf -std=c++17
# echo "done!"
# EOF
# """)
#         _, stdout, stderr = client.exec_command(f"chmod 770 {compile_script}")
#         # print(stdout.read(),stderr.read())
#         _, stdout, stderr = client.exec_command(f"{compile_script}")
#         print(stdout.read(),stderr.read())
