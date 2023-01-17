import pandas as pd
import numpy as np
import os, sys

root_dir = os.path.dirname(os.path.abspath(''))
if not root_dir in sys.path: sys.path.append(root_dir)

from empirical.readData_mat import *
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
            self.artificial_data(filePath)

        pass

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

    def xlsx_data(self, filePath='../../data/BuscheLab/2P_data.xlsx',sheets=None,population_keys=['mouse type']):

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

        # create dictionary to preformat data to be written into dataframe
        rate_dict = {}
        for sheet in sheets:        # iterate through specified sheets

            input_df = data.parse(sheet,header=list(range(len(population_keys)-1)))

            for col in input_df.columns:    # iterate through each column in sheet

                # create column keys from header-rows
                key = (sheet,*col) if isinstance(col,tuple) else (sheet,col)
                rate_dict[key] = pd.Series(input_df[col])

        self.rates = pd.DataFrame(rate_dict)
        column_names = [[str(entry) for entry in col] for col in self.rates.columns]
        self.rates.columns = pd.MultiIndex.from_tuples(column_names, names=population_keys)


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

        expand_df(population_names)

        # after all is done, sort dataframe such that indices appear at proper positions
        return data.sort_index(axis=1)


    def get_data_shape(self):

        """
            obtain maximum number of sublevels for each column level to create
            accordingly shaped tensor
        """

        population_names = list(self.rates.columns.names)

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



    def artificical_data(self,):

        self._num_layers = 1
        self._num_clusters = 2

        arr_rateWnt = [1.,2.,5.]
        arr_alpha_0 = [0.0,0.02,0.06]

        for L in range(arr_rateWnt):
            res = create_measures(L=1,S=[1,2],N=100,rerun=True,rateWnt=1.,alpha_0=0.02)
