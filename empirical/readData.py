import pandas as pd
import numpy as np
import os

from .readData_mat import *
from .readData_xls import *

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

    def mat_data(self, filePath='../../data/BuscheLab/spiking_data_for_modeling_with_depth.mat',plot=False):

        data = read_data(filePath=filePath,plot=plot)
        self._num_animals = len(data)
        self._num_layers = 3
        self._num_clusters = 2

        #types = [0 if sd["classification"]["spikes_genotype"][1]==114 else 1 for sd in data]

        self.modelShape = (self._num_animals,self._num_layers,self._num_clusters)
        self.types = ['WT','cTKO']
        self.type = np.array([0 if sd["classification"]["spikes_genotype"][1]==114 else 1 for sd in data])

        self.spikes_raw = np.zeros(self.modelShape + (0,))

        spikes_shape = self.spikes_raw.shape
        for animal_idx,sd in enumerate(data):

            for layer_idx in range(self._num_layers):

                for cluster_idx in range(self._num_clusters):

                    idx = (np.array(sd['layer'])==layer_idx) & (sd['cluster_idx']-1==cluster_idx)
                    nSpikes = idx.sum()
                    if nSpikes > spikes_shape[-1]:
                        self.spikes_raw = np.pad(self.spikes_raw,((0,0),(0,0),(0,0),(0,nSpikes-spikes_shape[-1])),'constant',constant_values=np.NaN)
                        spikes_shape = self.spikes_raw.shape

                    self.spikes_raw[animal_idx,layer_idx,cluster_idx,:nSpikes] = sd['rate'][idx,0]

        self.nMax = spikes_shape[-1]


        self.spikes = np.transpose(self.spikes_raw,(3,0,1,2))
        self.mask = ~np.isnan(self.spikes)
        self.spikes = self.spikes[self.mask]

    def xlsx_data(self, filePath='../../data/BuscheLab/2P_data.xlsx',sheets=None,include_silent=True):

        self._num_layers = 1
        self._num_clusters = 1

        data = pd.ExcelFile(filePath)
        sheets = sheets if sheets else data.sheet_names

        self._num_animals = 0
        spikes_shape = 0
        self.modelShape = (self._num_animals,self._num_layers,self._num_clusters)
        self.spikes_raw = np.zeros(self.modelShape + (spikes_shape,))
        self.types = sheets
        self.type = []

        for i,sheet in enumerate(sheets):

            read_spikes = pd.read_excel(filePath,sheet_name=sheet)
            read_spikes = read_spikes.to_numpy()

            nSpikes = read_spikes.shape[0]
            if nSpikes > spikes_shape:
                self.spikes_raw = np.pad(self.spikes_raw,((0,0),(0,0),(0,0),(0,nSpikes-spikes_shape)),'constant',constant_values=np.NaN)
                spikes_shape = nSpikes


            self.spikes_raw = np.pad(self.spikes_raw,((0,read_spikes.shape[1]),(0,0),(0,0),(0,0)),'constant',constant_values=np.NaN)

            for spikes in read_spikes.T:
                self.spikes_raw[self._num_animals,0,0,:nSpikes] = spikes
                self.type.append(i)
                self._num_animals += 1

        self.type = np.array(self.type)

        self.modelShape = (self._num_animals,self._num_layers,self._num_clusters)
        self.spikes = np.transpose(self.spikes_raw,(3,0,1,2))
        self.nMax = self.spikes.shape[0]
        if include_silent:
            self.spikes[self.spikes==0] = 1/600.

        self.mask = np.logical_and(~np.isnan(self.spikes),self.spikes>0)
        self.spikes = self.spikes[self.mask]

    def artificical_data(self,):

        self._num_layers = 1
        self._num_clusters = 2

        arr_rateWnt = [1.,2.,5.]
        arr_alpha_0 = [0.0,0.02,0.06]

        for L in range(arr_rateWnt):
            res = create_measures(L=1,S=[1,2],N=100,rerun=True,rateWnt=1.,alpha_0=0.02)
