import pandas as pd
import numpy as np
import os

from matplotlib import pyplot as plt

# root_dir = os.path.dirname(os.path.abspath(''))
# if not root_dir in sys.path: sys.path.append(root_dir)

# from .readData_mat import *
from darkMatter import darkMatter
from DM_theory.functions import *

#from empirical.readData_xls import *

class ModelParams:

    def __init__(self,mode,**kwargs):

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

        #types = [0 if sd["classification"]["spikes_genotype"][1]==114 else 1 for sd in data]

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


    def artificial_data(self,gamma,delta,nu_max,population_keys=['rates','mouse ID'],T=600.,N=100):
        
        self._num_layers = 1
        self._num_clusters = 2

        self.gamma = gamma
        self.delta = delta
        self.nu_max = nu_max

        self.N = N
        self.T = float(T)

        # arr_rateWnt = [1.]
        # arr_alpha_0 = [0.0,0.02,0.06]

        # self.spike_counts = pd.DataFrame(columns=arr_rateWnt)
        # print(self.rates)

        self.spike_counts = pd.DataFrame(
            columns=pd.MultiIndex(levels=[[],[]],codes=[[],[]],names=population_keys)
        )

        # for r,rate in enumerate(arr_rateWnt):
        options = {
            'L': 1,
            'P': 1,
            'S': [1],

            'tau_m': 0.01,
            'J_0': -1,
            'computation': {
                'N': self.N,
                'T': self.T,
                'draw_from_theory': 5,
                'draw_finite_time': 1,
            }
        }

        options['rateWnt'] = get_nu_bar(gamma=gamma,delta=delta,nu_max=nu_max)
        options['tau_I'] = [get_tau_I(nu_max=nu_max,tau_m=options['tau_m'])]
        options['alpha_0'] = get_alpha_0(gamma=gamma,delta=delta,nu_max=nu_max,tau_m=options['tau_m'],J_0=options['J_0'])

        print(f"inferred parameters: {options['rateWnt']=}, {options['tau_I']=}, {options['alpha_0']=}")

        res = darkMatter(steps=200,mode=1,options=options,cleanup=False,rerun=True,compile=False,logging=2)

        for d,rate_draw in enumerate(res['rates_T'][0,...,0].T):
            self.spike_counts = add_column_to_dataframe(self.spike_counts,rate_draw,d,population=f'{gamma=}',population_keys=population_keys)
        
        self.plot_rates()
        
        # return res
            # res = create_measures(L=1,S=[1,2],N=100,rerun=True,rateWnt=1.,alpha_0=0.02)


    def plot_rates(self):


        fig,ax = plt.subplots(1,1,figsize=(10,5))

        NU = np.linspace(0,self.nu_max,51)
        ## plot histogram of empirical or artificial rates
        ax.hist(self.spike_counts/self.T,bins=NU,density=True)



        ## plot underlying original distribution
        NU = np.linspace(0,self.nu_max,1001)
        ax.plot(NU,p_nu(NU,self.gamma,self.delta,self.nu_max),label='original distribution')

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

        #expand_df(population_names)

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
        return pd.concat([df,new_data_df],axis=1)