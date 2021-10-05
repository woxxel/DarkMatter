import h5py
from scipy.io import netcdf
import numpy as np

def mat2ncdf(fileName):
    
    # load data
    #fileMat = sio.loadmat(file_name)
    datFile = h5py.File(fileName)
    
    # write into dictionary
    data = {}
    
    data['N'] = datFile['data']['events'].shape[0]
    data['train'] = [[]]*data['N']
    data['N_AP'] = np.zeros(data['N'])
    data['rates'] = np.zeros(data['N'])
    
    freq = datFile[datFile['data']['frame_rate'][0,0]][0,0]
    
    data['T'] = datFile[datFile['data']['signal'][0,0]].shape[0]/freq
    print freq
    print datFile[datFile['data']['signal'][0,0]].shape
    print data['T']
    
    for n in range(data['N']):
        data['train'][n] = datFile[datFile['data']['events'][n,0]][0]/freq
        
        if type(data['train'][n]) == np.ndarray:
            data['N_AP'][n] = len(data['train'][n])
        else:
            data['N_AP'][n] = 0
        
        data['rates'][n] = data['N_AP'][n]/data['T']
    
    return data