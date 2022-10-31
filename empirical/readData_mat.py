import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import h5py

def compDiff(a,b,classify):
    comp = True
    for c in classify:
        if type(a[c])==np.ndarray:
            comp = comp and np.array_equal(a[c],b[c])
        else:
            comp = comp and a[c]==b[c]

    return comp


def plot_animal(spike_data,animal,ax,j=0):

    i=0
    for sd in spike_data:
        if sd['classification']['spikes_animal_idx'] == animal:

            for n in range(len(sd['neuron_idx'])):
                nAP = len(sd['spike_times'][n])
                col = 'gray'
                if sd['cluster_idx'][n] == 1.:
                    col = 'k'
                elif sd['cluster_idx'][n] == 2.:
                    col = 'r'
                ax[i][0].plot(sd['spike_times'][n],np.ones(nAP)*n,'.',markersize=2,color=col)
            ax[i][0].set_title(f'animal {int(animal)}' \
                                + f'{"l" if sd["classification"]["spikes_hemisphere_idx"]==1 else "r"}' \
                                + f' {"cre" if sd["classification"]["spikes_genotype"][1]==114 else "cTKO"}')
            ax[i][1].hist(sd['rate'][sd['cluster_idx']==1,j],np.linspace(0,5,51),facecolor='k',alpha=0.8)
            ax[i][1].hist(sd['rate'][sd['cluster_idx']==0,j],np.linspace(0,5,51),facecolor='gray',alpha=0.8)
            ax[i][1].hist(sd['rate'][sd['cluster_idx']==2,j],np.linspace(0,5,51),facecolor='r',alpha=0.8)

            i+=1
    plt.tight_layout()

def read_data(filePath='../data/BuscheLab/spiking_data_for_modeling_with_depth.mat', classify = ['spikes_animal_idx','spikes_hemisphere_idx','spikes_genotype'],plot=False):

    """

        filePath: string
            path to empirical data
        classify: list(strings)
            list, specifying data along which to cluster
        plt: bool [False]
            specifies, whether to plot data

    """

    ### read data from file
    f = h5py.File(filePath,'r')
    data = {}
    for key in f.keys():
        data[key] = np.array(f.get(key))

    ### preprocess data where needed
    data['spikes_genotype'] = []
    for spike_genotype_ref in data['spikes_genotype_idx']:
        data['spikes_genotype'].append(np.array(f[spike_genotype_ref[0]]))
    data['spikes_animal_idx'] = data['spikes_animal_idx'].flatten()
    data['spikes_cluster_idx'] = data['spikes_cluster_idx'].flatten()
    data['spikes_hemisphere_idx'] = data['spikes_hemisphere_idx'].flatten()
    data['spikes_depth_rel2L5'] = data['spikes_depth_rel2L5'].flatten()

    ### Compart neurons into layers
    layer_border = [-1000,-250,200,1000]
    depth = np.array(data['spikes_depth_rel2L5'])
    data['layer'] = np.zeros((len(depth)),'int')
    for l in range(len(layer_border)-1):
        #print((layer_border[l-1] < depth) & (depth < layer_border[l]))
        data['layer'][(layer_border[l] <= depth) & (depth < layer_border[l+1])] = l
        #print(l,layer_border[l],layer_border[l+1], sum(data['layer']==l))

    ### clustering data
    spike_data = []
    for n in range(len(data['spikes_animal_idx'])):

        ### create classification dictionary
        tmp_class = {}
        for c in classify:
            tmp_class[c] = data[c][n]

        ### search, whether this set of data already exists
        idx = None
        for i,sd in enumerate(spike_data):
            if compDiff(tmp_class,sd['classification'],classify):
                idx = i
                break

        ### create new cluster if idx not found
        if idx==None:
            idx = len(spike_data)
            spike_data.append({
                'classification': tmp_class,
                'spike_times': [],
                'neuron_idx': [],
                'cluster_idx': [],
                'layer': [],
                'depth': [],
            })

        spike_data[idx]['spike_times'].append(np.array(f[data['spikes_timestamps'][n][0]])[0])
        spike_data[idx]['neuron_idx'].append(n)
        spike_data[idx]['cluster_idx'].append(data['spikes_cluster_idx'][n])
        spike_data[idx]['layer'].append(data['layer'][n])
        spike_data[idx]['depth'].append(data['spikes_depth_rel2L5'][n])


    ### some processing

    ### sort neurons according to type
    sort_keys = ['neuron_idx','cluster_idx','spike_times']
    for sd in spike_data:
        sort_idx = np.argsort(sd['cluster_idx'])

        sd_tmp = {}
        for key in sort_keys:
            data_tmp = []
            for n in sort_idx:
                data_tmp.append(sd[key][n])
            # if not (key=='spike_times'):
            #     print('data')
            #     print(data_tmp)
            sd[key] = data_tmp
        sd['cluster_idx'] = np.array(sd['cluster_idx'])

    ### calculate some statistics
    ## firing rate
    for sd in spike_data:
        sd['N'] = len(sd['neuron_idx']) # total number of recorded neurons
        sd['g'] = np.array([np.sum(np.array(sd['cluster_idx'])==a) for a in range(3)])/sd['N'] # fraction of neuron belonging to cluster a

        sd['T'] = max([max(s) for s in sd['spike_times']]) # recording time (estimated from very last spike encountered in this animal)

        nIntervals = 6
        dT = sd['T']/nIntervals
        sd['rate'] = np.zeros((sd['N'],nIntervals+1))
        for n,s in enumerate(sd['spike_times']):
            sd['rate'][n,0] = len(s)/sd['T']
            for i in range(nIntervals):
                sd['rate'][n,i+1] = np.sum((s>i*dT) & (s<(i+1)*dT))/dT

    if plot:
        nAnimals = len(np.unique(data['spikes_animal_idx']))

        r = np.zeros((nAnimals,2))
        N = np.zeros((nAnimals,2))

        for sd in spike_data:
            a = int(sd['classification']['spikes_animal_idx'])-1
            h = int(sd['classification']['spikes_hemisphere_idx'])-1
            r[a,h] = sd['g'][2]/np.sum(sd['g'][1:])
            N[a,h] = sd['N']

        fig, ax = plt.subplots(2,2,figsize=(10,8))
        ax[0][0].plot(N[:,0],r[:,0],'ro',label='left')
        ax[0][0].plot(N[:,1],r[:,1],'ko',label='right')
        ax[0][0].plot(N.T,r.T,'k-')

        #print('fraction of inhibitory neurons: ', r)
        plt.setp(ax[0][0],ylabel='fraction of inhibitory neurons',xlabel='total number of neurons')
        ax[0][0].legend()


        N = np.zeros((nAnimals,3,3))

        for a in range(nAnimals):
            for h in range(2):
                for l in range(3):
                    for c in range(3):
                        for sd in spike_data:
                            if int(sd['classification']['spikes_animal_idx'])-1==a:
                                idx = (np.array(sd['layer'])==l) & (sd['cluster_idx']==c)
                                N[a,l,c] += idx.sum()

        ax[0][1].plot(N[:,0,1],N[:,0,2],'o',color=[0,0,0.2],label='L1')
        ax[0][1].plot(N[:,1,1],N[:,1,2],'o',color=[0,0,1],label='L2')
        ax[0][1].plot(N[:,2,1],N[:,2,2],'o',color=[0.7,0.7,1],label='L3')
        plt.setp(ax[0][1],xlabel='excitatory',ylabel='inhibitory')

        ax[0][1].legend()


        _ = ax[1][0].hist(data['spikes_depth_rel2L5'],np.linspace(-1000,1000,101))
        plt.setp(ax[1][0],xlabel='depth',ylabel='# neurons')
        plt.show(block=False)

        fig,ax = plt.subplots(4,5,figsize=(15,8))
        for i,animal in enumerate(spike_data):
            _ = ax[i//5][i%5].hist(animal['depth'],np.linspace(-1000,1000,51))


        plt.show(block=False)

        animal = 2
        fig,ax = plt.subplots(2,2,figsize=(12,4))

        def change_slider(animal):
            # print(animal)
            plot_animal(spike_data,animal,ax)

        axamp = plt.axes([0.02, .6, 0.4, 0.02])
        plot_animal(spike_data,animal,ax,4)
        slider = Slider(axamp, 'animal', 1, 9, valinit=animal,orientation='horizontal')
        slider.on_changed(change_slider)
        plt.show(block=False)

    return spike_data
