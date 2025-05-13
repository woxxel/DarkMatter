# import pandas as pd
import numpy as np
import os,sys

from matplotlib import pyplot as plt

root_dir = os.path.dirname(os.path.abspath(''))
if not root_dir in sys.path: sys.path.append(root_dir)

# from .readData_mat import *
# from darkMatter import darkMatter
from DM_theory.functions import get_nu_bar,get_tau_I,get_alpha_0

from general_utils.parameters import set_options
from inference.utils.utils import *

# from empirical.readData_xls import *

class ModelParams:

    def __init__(self,mode,path_results=None,suffix='',**kwargs):

        if suffix !='':
            self.suffix = f"_{suffix}" if not suffix.startswith('_') else suffix
        else:
            self.suffix = ''

        self.path_results = path_results
        if mode=='empirical':

            filePath = kwargs['filePath']
            ext = os.path.splitext(filePath)[1]
            if ext=='.mat':
                self.mat_data(filePath)
            elif ext=='.xlsx':
                self.xlsx_data(**kwargs)

        if mode=='artificial':
            self.artificial_data(**kwargs)

        # pass

    # def initialize_rate_df(self,mouse_types,nMice,n_layer=1,n_cluster=1):
    #
    #     """ initializing an empty pandas-dataframe to store neuron firing rate data
    #
    #         Input parameters:
    #             mouse_types: list(str)
    #                 list of mouse types
    #             nMice:       list(int)
    #                 number of mice per mouse type
    #             n_layer:     int [default: 1]
    #                 number of layers in the data
    #             n_cluster:   int [default: 1]
    #                 number of neuron cluster types (exc., inh., ...)
    #
    #         Output:
    #             df: pandas Dataframe
    #                 empty dataframe with multicolumn structure for storing
    #                 firing rates
    #     """
    #
    #     ## checking if input makes sense
    #     assert isinstance(mouse_types,list) and isinstance(nMice,list), f'type of mouse_types ({type(mouse_types)}) and nMice ({type(nMice)}) required to be lists.'
    #     assert len(mouse_types)==len(nMice), f'number of mouse types ({len(mouse_types)}) and mouse numbers ({len(nMice)}) need to agree.'
    #
    #     ## creating the column labels
    #     col_names = []
    #     for type_idx,mouse_type in enumerate(mouse_types):
    #         for mouse_ID in range(nMice[type_idx]):
    #             for layer in range(n_layer):
    #                 for cluster in range(n_cluster):
    #                     col_names.append((mouse_type,mouse_ID,layer,cluster))
    #
    #     ## creating and returning the empty dataframe
    #     df = pd.DataFrame(columns=col_names)
    #     df.columns = pd.MultiIndex.from_tuples(df.columns, names=['mouse type','mouse ID','layer','cluster'])
    #     return df

    def mat_data(self, filePath='../../data/BuscheLab/spiking_data_for_modeling_with_depth.mat',plot=False):

        data = read_data(filePath=filePath,plot=plot)
        self._num_animals = len(data)
        self._num_layers = 3
        self._num_clusters = 2

        # types = [0 if sd["classification"]["spikes_genotype"][1]==114 else 1 for sd in data]

        self.modelShape = (self._num_animals,self._num_layers,self._num_clusters)
        self.types = ['WT','cTKO']
        self.type = np.array([0 if sd["classification"]["spikes_genotype"][1]==114 else 1 for sd in data])

        self.rates_raw = np.zeros(self.modelShape + (0,))
        self.nMice = np.zeros(len(self.types))

        rates_shape = self.rates_raw.shape
        for type_idx, animal_type in enumerate(self.types):

            animal_idx = 0

            for sd in data[self.type==type_idx]:

                for layer_idx in range(self._num_layers):

                    for cluster_idx in range(self._num_clusters):

                        idx = (np.array(sd['layer'])==layer_idx) & (sd['cluster_idx']-1==cluster_idx)
                        nRates = idx.sum()
                        if nRates > rates_shape[-1]:
                            self.rates_raw = np.pad(self.rates_raw,((0,0),(0,0),(0,0),(0,nRates-rates_shape[-1])),'constant',constant_values=np.NaN)
                            rates_shape = self.rates_raw.shape

                        self.rates_raw[type_idx,animal_idx,layer_idx,cluster_idx,:nRates] = sd['rate'][idx,0]
                animal_idx += 1

            self.nMice[type_idx] = animal_idx

        self.nMax = rates_shape[-1]

        self.rates = np.transpose(self.rates_raw,(3,0,1,2))
        self.mask = ~np.isnan(self.rates)
        self.rates = self.rates[self.mask]

    def xlsx_data(self, filePath='../../data/BuscheLab/2P_data.xlsx',sheets=None,population_keys=['mouse type','mouse_ID']):

        """
            function to turn xlsx sheets into a properly structured pandas-DataFrame
            and write results to self.rates

            Input:
                filePath: string
                    (relative/absolute) path to the file
                sheets: list(string)
                    list of sheets to be read from the file. None (default) reads all
                population_keys: list(string)
                    Meta-tags to determine names of population. First is descriptor
                    of sheetnames followed by descriptor for title-rows of excel.
                    len(population_keys)-1 equals the number of headlines expected
                    in file and len(population)_keys) is the number of multicolumn
                    indices in resulting pandas-DataFrame
        """
        assert len(population_keys)>=1, 'at least one key for population distinction has to be provided'

        data = pd.ExcelFile(filePath)
        sheets = sheets if sheets else data.sheet_names

        self.rates = pd.DataFrame(
            columns=pd.MultiIndex(levels=[[],[]],codes=[[],[]],names=population_keys)
        )

        for sheet in sheets:        # iterate through specified sheets

            input_df = data.parse(sheet,header=list(range(len(population_keys)-1)))
            for col in input_df.columns:    # iterate through each column in sheet

                self.rates = add_column_to_dataframe(self.rates,new_data=np.array(input_df[col]),name=col,population=sheet,population_keys=population_keys)

        self.get_data_shape()

        self.types = np.unique(self.rates.columns.get_level_values(population_keys[0]))

        self.animal_mask = np.zeros(np.prod(self.data_shape),dtype='bool')
        for i,c in enumerate(self.types):
            self.animal_mask[i*self.data_shape[1]:i*self.data_shape[1]+self.rates[c].shape[1]] = True

        self.T = 1200

    def artificial_data(self,parameter,population_keys=['rates','mouse ID'],T=600.,N=100,nAnimals=5,plot=False):

        self._num_layers = 1
        self._num_clusters = 2

        self.params = parameter
        self.two_pop = 'p' in self.params

        self.N = N
        self.T = float(T)

        # arr_rateWnt = [1.]
        # arr_alpha_0 = [0.0,0.02,0.06]

        # self.spike_counts = pd.DataFrame(columns=arr_rateWnt)
        # print(self.rates)

        # self.spike_counts = pd.DataFrame(
        #     columns=pd.MultiIndex(levels=[[],[]],codes=[[],[]],names=population_keys)
        # )

        options = set_options(nI=2 if self.two_pop else 1,nE=0)
        options['eps'] = 0.

        # for r,rate in enumerate(arr_rateWnt):
        options['computation'] = {
                'N': self.N,
                'T': self.T,
                'draw_from_theory': nAnimals,
                'draw_finite_time': 1,
            }

        options['rateWnt'] = []
        options['tau_I'] = []
        options['alpha_0'] = []

        for m in range(len(self.params["distr"])):
            options["rateWnt"].append(
                get_nu_bar(
                    gamma=self.params["distr"][m]["gamma"],
                    delta=self.params["distr"][m]["delta"],
                    nu_max=self.params["distr"][m]["nu_max"],
                )
            )
            options["tau_I"].append(
                get_tau_I(
                    nu_max=self.params["distr"][m]["nu_max"], tau_m=options["tau_M"]
                )
            )
            options["alpha_0"].append(
                get_alpha_0(
                    gamma=self.params["distr"][m]["gamma"],
                    delta=self.params["distr"][m]["delta"],
                    nu_max=self.params["distr"][m]["nu_max"],
                    tau_m=options["tau_M"],
                    J_0=options["J0"][0],
                )
            )
        # options['mode_selfcon'] = 0
        # print(options)
        # print('\n')
        string_input_params = "input parameters: "
        for m in range(len(self.params["distr"])):
            for key in ['gamma','delta','nu_max']:
                string_input_params += f"{self.params['distr'][m][key]=}, "
        print(string_input_params)

        print(f"inferred parameters: {options['rateWnt']=}, {options['tau_I']=}, {options['alpha_0']=}")

        if np.isnan(options['rateWnt']).any() or np.isnan(options['tau_I']).any() or np.isnan(options['alpha_0']).any():
            print('inferred parameters contain NaNs - breaking!')
            self.rates = False
            return

        self.spike_counts = np.full((int(self.N), nAnimals), np.nan)
        for a in range(nAnimals):
            samples, samples_T = draw_samples(f_target=p_nu,params=self.params,n_samples=self.N,T=self.T,tolerance=0.001,plot=False,save=False)

            self.spike_counts[:,a] = samples_T
            # self.spike_counts = add_column_to_dataframe(self.spike_counts,samples_T,a,population=f"{self.params['gamma_1']=}",population_keys=population_keys)

        # res = darkMatter(steps=200,mode=1,options=options,path=self.path_results,cleanup=True,rerun=True,compile=True,logging=2,suffix=self.suffix)
        # return res
        # print(res)

        # print(f'{two_pop=}')

        # rate_lims = [-4,1.5]

        # self.res = res
        # print(res['rates_T'].shape)
        # fig,ax = plt.subplots(2,5,figsize=(10,15))
        # for d,_rate_draw in enumerate(np.transpose(res['rates_T'][...,0],(2,0,1))):
        #     # print(_rate_draw.shape)
        #     # should take rate_T[0 and 1], with fraction according to 'p'
        #     if two_pop:
        #         rate_draw = np.concatenate([
        #             np.random.choice(_rate_draw[0,:],int(self.params['p']*_rate_draw.shape[1])),
        #             np.random.choice(_rate_draw[1,:],int((1-self.params['p'])*_rate_draw.shape[1]))
        #         ],axis=0)
        #         # print(_rate_draw[0,:]/self.T)
        #         # ax[0][d].hist(_rate_draw[0,:]/self.T,bins=np.logspace(rate_lims[0],rate_lims[1],51),density=True,color='tab:blue',alpha=0.6)
        #         # ax[0][d].hist(_rate_draw[1,:]/self.T,bins=np.logspace(rate_lims[0],rate_lims[1],51),density=True,color='tab:red',alpha=0.6)
        #         # ax[0][d].set_yscale('log')
        #         # ax[0][d].set_xscale('log')
        #         # plt.setp(ax[0][d],ylim=[10**-4,10**1],xlim=[10**rate_lims[0],10**rate_lims[1]])
        #         # # print(NU)
        #         # NU = np.linspace(0,self.params['nu_max'],10**6+1)
        #         # p_nu_1 = p_nu(NU,self.params,two_pop=False)
        #         # p_nu_2 = p_nu(NU,{'gamma':self.params['gamma2'],'delta':self.params['delta2'],'nu_max':self.params['nu_max2']},two_pop=False)
        #         # ax[0][d].plot(NU,p_nu_1,color='tab:blue',linestyle='--')
        #         # ax[0][d].plot(NU,p_nu_2,color='tab:red',linestyle='--')

        #         # ax[1][d].hist(rate_draw/self.T,bins=np.logspace(rate_lims[0],rate_lims[1],51),density=True,color='tab:green',alpha=0.6)
        #         # ax[1][d].plot(NU,p_nu_1*self.params['p']+p_nu_2*(1-self.params['p']),color='tab:green',linestyle='--')
        #         # ax[1][d].set_yscale('log')
        #         # ax[1][d].set_xscale('log')
        #         # plt.setp(ax[1][d],ylim=[10**-4,10**1],xlim=[10**rate_lims[0],0.5*10**rate_lims[1]])

        #     else:
        #         rate_draw = _rate_draw[0,:]

        #     # print(rate_draw,rate_draw.shape)
        #     self.spike_counts = add_column_to_dataframe(self.spike_counts,rate_draw,d,population=f"{self.params['gamma']=}",population_keys=population_keys)

        # plt.show(block=False)
        self.rates = self.spike_counts/self.T
        # print(self.rates)
        if plot:
            self.plot_rates(key=f"{self.params['distr'][0]['gamma']=}")

        # return res
        # res = create_measures(L=1,S=[1,2],N=100,rerun=True,rateWnt=1.,alpha_0=0.02)

    def plot_rates(self,key=None,param_in=None):

        fig,ax = plt.subplots(1,2,figsize=(10,5))

        # gamma = gamma if gamma else self.params['gamma']
        # delta = delta if delta else self.params['delta']
        # nu_max = nu_max if nu_max else self.params['nu_max']
        # parameters = [self.params]
        # if param_in:
        # parameters.append(param_in)
        param = param_in if param_in else self.params

        # rates = self.rates[key] if key else self.rates
        rates = self.rates

        xlim = 25.

        print(f"{self.two_pop=}")

        bins = np.linspace(0, param["distr"][0]["nu_max"], 101)

        NU_log = np.logspace(-20, np.log(param["distr"][0]["nu_max"]), 10001)
        p_NU = p_nu(NU_log,param)
        p_NU_cum = np.nancumsum(p_nu(NU_log,param)[:-1]*np.diff(NU_log))

        ## plot histogram of empirical or artificial rates
        ax[0].hist(rates,bins=bins,density=True)
        # ax[0].set_xscale('log')

        ## plot underlying original distribution
        ax[1].hist(rates,bins=bins,density=True,cumulative=True,histtype='step')
        # print(parameters)
        colors = ['k','r']
        # NU = np.linspace(0,xlim,10**8+1)
        # for i,param in enumerate(parameters):
        # print(p_NU)
        ax[0].plot(NU_log,p_NU,label='original distribution')

        # bins = 10**np.linspace(-4,2,101)
        # p_NU[-1] = 0
        # p_NU_cum = np.nancumsum(p_NU)
        # p_NU_cum /= np.nanmax(p_NU_cum)
        # print(p_NU,p_NU_cum)

        ax[1].plot(NU_log[:-1],p_NU_cum,label='original distribution',color='k',linestyle='--')

        plt.setp(ax[0],xlim=[10**(-4),xlim],ylim=[0,2])
        plt.setp(ax[1],xlim=[10**(-4),xlim],ylim=[0,1.1])
        plt.show(block=False)

    def regularize_rates(self):
        """
            this method is used to create a regular dataframe with the same number
            of sub-levels on each level, to allow proper interaction
        """

        population_names = list(self.rates.columns.names)
        max_shape = self.get_data_shape()

        data = self.rates.copy()

        def expand_df(levels,selector=[]):

            """
                manipulates self.rates DataFrame to have regular shape

                Input:
                    levels: list(string)
                        remaining levels to dive into
                    selector: list(string)
                        keys to navigate to current section of data
            """

            # obtain lower level keys
            this_level = levels[0]
            if len(selector):
                keys = np.unique(list(data[selector].columns.get_level_values(this_level)))
            else:
                keys = np.unique(list(data.columns.get_level_values(this_level)))

            # if some data is missing, create dummy entries
            # if len(keys) < max_shape[this_level]:
            for i in range(max_shape[this_level] - len(keys)):
                print(f'expanding @ lvl {this_level} with keys: {keys}, selectors: {selector}')

                name = f'dummy_{i}'
                j=1
                while name in keys:
                    name = f'dummy_{i+j}'
                    j += 1

                key = (*selector,*[name]*len(levels))
                data[key] = np.NaN

                keys = np.append(keys,name)

            if len(levels)<=1:
                return

            # dive deeper
            for key in keys:
                expand_df(levels[1:],(*selector,key))

        # expand_df(population_names)

        # after all is done, sort dataframe such that indices appear at proper positions
        return data.sort_index(axis=1)

    def get_data_shape(self):

        """
            obtain maximum number of sublevels for each column level to create
            accordingly shaped tensor
        """

        population_names = list(self.rates.columns.names)
        print(population_names)

        max_shape = {}
        for name in population_names:
            max_shape[name] = 0

        def get_max_dim(data,levels):

            this_level = levels[0]
            keys = np.unique(list(data.columns.get_level_values(this_level)))

            max_shape[this_level] = max(max_shape[this_level],len(keys))

            if len(levels)<=1:
                return

            for key in keys:
                get_max_dim(data[key],levels[1:])

        get_max_dim(self.rates,population_names)
        # max_shape['val1'] = 2

        self.data_shape = list(max_shape.values())

        return max_shape


def add_column_to_dataframe(df,new_data,name,population=None,population_keys=None):

    """
        function to add a column to a pandas dataframe
        with multicolumn structure

        Input:
            df: pandas DataFrame
                DataFrame to which the column is added
            data: np.array
                data to be added
            name: string
                name of the column
            population: string
                name of the population (if provided, introduces hierarchical level)
    """

    # construct index to multicolumn
    if population:
        idx = (population,name)
    else:
        idx = name

    new_data_df = pd.DataFrame(
        new_data,
        columns=pd.MultiIndex.from_tuples([idx],names=population_keys)
    )

    if (df is None) or df.empty:
        return new_data_df
    else:
        if len(new_data_df) > len(df):
            df = df.reindex(new_data_df.index)
        return pd.concat([df, new_data_df], axis=1)


def draw_samples(f_target,params,n_samples=1000,T=1200,tolerance=0.001,plot=False,save=False):

    ## define bounding function
    # f_envelope = lambda nu,M : M* 1./nu

    high = params["distr"][0]["nu_max"]
    ## finding the lower bound for the bounding function by requiring the CDF(g(nu),0,b) < 0.001
    i = -1
    while True:
        low = 10**i

        res,err = quad(f_target,0,low,args=(params),points=np.logspace(-3,0,4))
        if res < tolerance:
            break

        i -= 1
    # print(f'lower bound: {low=}')

    def normalized_envelope_density(x, low, high):
        return 1 / (x * np.log(high / low))

    def sample_from_normalized_envelope(low, high):
        u = np.random.uniform(0, 1)
        return low * (high / low)**u  # Sample from g(x) normalized over [b, a]

    def rejection_sampling(n_samples, M, low, high):
        samples = []
        while len(samples) < n_samples:
            nu = sample_from_normalized_envelope(low, high)
            v = np.random.uniform(0, 1)
            if v <= f_target(nu,params) / (M * normalized_envelope_density(nu, low, high)):
                samples.append(nu)
        return np.array(samples)

    # define M such that M*g(nu) > f(nu) for all nu
    nu = np.logspace(np.log10(low),np.log10(high),10**3+1)
    M = np.ceil(np.nanmax(f_target(nu,params) / normalized_envelope_density(nu,low,high)))

    samples = rejection_sampling(n_samples, M, low, high)
    samples_T = np.random.poisson(samples*T,samples.shape)

    if plot:
        fig,ax = plt.subplots(1,2,figsize=(6,3))
        ax[0].hist(samples,bins=np.logspace(np.log10(low),np.log10(high),101),density = True,label='samples')
        ax[0].plot(nu, f_target(nu,params),'k-',label='target f')
        ax[0].plot(nu, M*normalized_envelope_density(nu,low,high),'r-',label='envelope g')
        ax[0].axvline(low,color='r',linewidth=0.5,linestyle='--',label='lower bound')
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].legend(bbox_to_anchor=(1.05, 1.05), loc='upper right')

        ax[1].hist(samples,bins=np.linspace(0,high,101),density = True)
        ax[1].plot(nu, f_target(nu,params),'k-')
        ax[1].plot(nu, M*normalized_envelope_density(nu,low,high),'r-')
        plt.setp(ax[1],xlim=[0,5],ylim=[0,1])

        for axx in ax:
            plt.setp(axx,xlabel='rate [Hz]',ylabel='density')
            axx.spines[['top','right']].set_visible(False)

        plt.tight_layout()
        if save:
            plt.savefig(f'./figures/rejection_sampling_example.png')
        plt.show(block=False)

    return samples, samples_T
