import numpy as np
import os
from netCDF4 import Dataset, stringtochar

"""
    This file contains the wrapping structure to call and read the c++ script,
    calling the computations for evaluating firing rate distributions and a number
    of parameters thereof

    steps:      int (100)
        number of steps taken along parameters that should be varied

    options:    dict ({})
        holds any number of parameters to override default values (see in set_model,
        set_simulation) with according (key,value) pairs. See specific description
        of parameters below

    rerun: bool (False)
        defines, whether existing result files with according parameters should
        be overriden (True), or used as output file (False)

    compile:    bool (False)
        defines, whether the c++ code should be compiled before running (only necessary,
        when changes to the code are made)
"""

def darkMatter(steps=100,options={},rerun=False,compile=False):

    # set model parameters
    fileModel = './data/modPara.nc'
    fileSim = './data/simPara.nc'

    model = set_model(fileModel,options)
    sim, sv_str = set_simulation(fileSim,options,steps)

    # fileResults = './data/results.nc'
    fileResults = './data/results_%s.nc' % (sv_str)

    if not os.path.exists(fileResults) or rerun:
        if (compile):
            compile_code()

        run_str = './theory/sharkfins ' + fileModel + ' ' + fileSim + ' ' + fileResults
        os.system(run_str)

    # obtain data from simulation
    ncid = Dataset(fileResults, "r")

    results = {}
    for key in ncid.variables:
        results[key] = ncid.variables[key][:]
    ncid.close()

    return results



class parameters:

    # class to support writing datafiles
    def __init__(self,name):
        # print('Building %s class'%name)
        self.paras = []

    def set_para(self,key,default,option=None,register=True,sz=None):

        if option:
            val = option[key] if key in option.keys() else default
        else:
            val = default
        if sz:
            if not (type(val)==list or type(val)==np.ndarray):
                val = np.full(sz,val)
                # val *= np.ones(sz)
            else:
                assert len(val)==sz, f'Please provide values for {key} of dimension 1 or {sz}'

        setattr(self,key,val)
        # if key in option.keys():
        # else:
        #     setattr(self,key,default)

        if register:
            self.paras.append(key)

    def prepare_sim_paras(self,steps):   # only works in conjunction with "set_simulation"
        for key in self.paras:
            val = getattr(self,key)
            assert type(val) == list, 'Please specify all iteration parameters as lists!'
            if len(val) == 2:
                val = np.linspace(val[0],val[-1],steps+1)[1:]#.astype('float')
            elif len(val) > 2:
                val = np.array(val)
            else:
                val = np.array(val).astype('float')
            setattr(self,key,val)

        nChar = len(max(self.paras,key=len))
        self.order = np.array(self.paras,dtype='S'+str(nChar))

    def print_parameter(self):

        for key in self.paras:
            print(key, getattr(self,key))

def set_model(fileModel,options):

    # defining model constants (all of them!)
    model = parameters('model')

    ## global parameters
    model.set_para('L',1,options)
    L = getattr(model,'L')
    assert not (L == len(options['S'])), f'Number of populations must be twice the number of layers ({L}), Npop={len(options["S"])} given'
    P = 2       # populations per layer
    nP = 2*L    # total number of populations

    model.set_para('S',[1,2],options)
    S = getattr(model,'S')
    layer_L_idx = []
    pop_L_idx = []
    pop_P_idx = []
    psp_L_idx = []
    psp_P_idx = []
    psp_S_idx = []

    for l in range(L):
        layer_L_idx.append(l)
        for p in range(P):
            pop_L_idx.append(l)
            pop_P_idx.append(p)
            for s in range(S[p+2*l]):
                psp_L_idx.append(l)
                psp_P_idx.append(p)
                psp_S_idx.append(s)

    nS = sum(S)
    # print('nS:',nS)
    # print('indices')
    # print(layer_L_idx)
    # print(pop_L_idx)
    # print(pop_P_idx)
    # print(psp_L_idx)
    # print(psp_P_idx)
    # print(psp_S_idx)


    ### parameters for layers
    model.set_para('layer_L_idx',layer_L_idx)

    model.set_para('eps',0.,options,sz=L)
    model.set_para('eta',0.,options,sz=L)
    model.set_para('J0_l',1.,options,sz=L)              # synaptic strength of inter-layer interactions
    model.set_para('kappa',1.,options,sz=L)            # connectivity ratio (should have dim=Npop)

    ### parameters for each population
    model.set_para('pop_L_idx',pop_L_idx)
    model.set_para('pop_P_idx',pop_P_idx)

    model.set_para('I_ext',1,options,sz=nP)           # external, constant drive
    model.set_para('rateWnt',1.,options,sz=nP)           # external, constant drive
    model.set_para('alpha_0',0.,options,sz=nP)           # external, constant drive
    model.set_para('tau_M',0.01,options,sz=nP)        # membrane timeconstant in sec
    model.set_para('tau_n',0.,options,sz=nP)
    model.set_para('J0',-1.,options,sz=nP)             # base synaptic strength

    model.set_para('drive',0,options,sz=L)            # external drive (1=on, 0=off)

    ### parameters for each synapse set
    model.set_para('psp_L_idx',psp_L_idx)
    model.set_para('psp_P_idx',psp_P_idx)
    model.set_para('psp_S_idx',psp_S_idx)

    model.set_para('tau_I',[0.005,0.2,0.03],options,sz=nS)
    model.set_para('tau_norm',1.,options,sz=nS)

    for tau_para in ['tau_I','tau_norm']:
        assert len(getattr(model,tau_para)) == nS, 'Make sure to provide the same number of parameters for each variable of synaptic timeconstants'
        assert len(getattr(model,tau_para)) >= nP, 'Number of defined synaptic timeconstants needs to be at least as large as the number of populations'

    # for p in range(Npop):
    #     if (model.drive[p] > 0):
    #         model.set_para('tau_0',1.,options,sz=Npop)         # external timeconstant
    #         model.set_para('J_0',1.,options,sz=Npop)           # external synaptic weight
    #         model.set_para('K_0',1,options,sz=Npop)            # external in-degree
    #     else:
    #         model.set_para('tau_0',0.,options,sz=Npop)         # external timeconstant
    #         model.set_para('J_0',0.,options,sz=Npop)           # external synaptic weight
    #         model.set_para('K_0',0,options,sz=Npop)            # external in-degree


    model.set_para('I_alpha',1.,options)
    model.set_para('I_beta',1.,options)

    ncid = Dataset(fileModel, "w");#, format="NETCDF4")
    ncid.createDimension('one',1)
    for key in model.paras:

        val = getattr(model,key)
        varType,varDim = prepare_ncid_var(ncid,model,key,val)
        Var = ncid.createVariable(key,varType,varDim)
        Var[:] = val
    ncid.close()

    # model.print_parameter()

    return model

def set_simulation(fileSim,options,steps):

    # defining simulation parameters (only the ones, which are being iterated)
    sim = parameters('simulation')

    for key in options['simulation']:
        if key in ['sim_prim','sim_sec']:
            sim.set_para(key,[],options['simulation'],register=False)
        else:
            sim.set_para(key,[],options['simulation'])

    # sim.set_para('order',sim.paras,{},register=False)

    # sim.set_para('rateWnt',[0.,20.],options['simulation'])
    # sim.set_para('alpha_0',[0.,0.16],options['simulation'])
    # sim.set_para('tau_G',[0.005],options)
    # sim.set_para('n',[0.],options)
    #
    # sim.set_para('I_alpha',[1],options)
    # sim.set_para('I_beta',[1],options)

    sim.set_para('mode_calc',0,options,register=False)        # 0=exact, 1=approx
    sim.set_para('mode_stats',0,options,register=False)       # 0=general stats, 1=...

    # order = ['rateWnt','alpha_0','tau_G','n','eps','eta','I_alpha','I_beta']
    # sim.set_para('order',order,options)
    sim.prepare_sim_paras(steps)

    sv_str = 'mode=%d_steps=%d' % (sim.mode_stats,steps)
    sv_str += '_iter'
    for key in sim.paras:
        val = getattr(sim,key)
        if len(val) > 1:
            sv_str += '_%s=%g' % (key,val[-1])

    sv_str += '_const'
    for key in sim.paras:
        val = getattr(sim,key)
        if (len(val) == 1 or val[0]==val[-1]):
            sv_str += '_%s=%g' % (key,val[0])

    #write simulation parameters to netcdf file
    ncid = Dataset(fileSim, "w")
    ncid.createDimension('one',1)

    sim.paras.extend(['mode_calc','mode_stats','sim_prim','sim_sec','order'])
    # print(sim.__dict__)
    for key in sim.paras:

        val = getattr(sim,key)
        varType,varDim = prepare_ncid_var(ncid,sim,key,val)
        if varType=='S':        ## somewhat complicated passing over readable strings to netcdf
            nChar = len(max(val,key=len))
            ncid.createDimension('charSz',nChar)
            Var = ncid.createVariable(key,'%s1'%(varType),varDim+('charSz',))
            Var._Encoding = 'ascii'
            Var[:] = stringtochar(val)
        else:
            Var = ncid.createVariable(key,varType,varDim)
            Var[:] = val

    ncid.close()

    # sim.print_parameter()

    return sim, sv_str

def get_type(val):
    if type(val)==np.int_ or type(val)==int:
        varType = np.dtype('int32').char #'i' #
    elif type(val)==np.bytes_ or type(val)==bytes:
        varType='S'
    else:
        varType = np.dtype('float64').char
    return varType

def prepare_ncid_var(ncid,sim,key,val):

    if type(val) == np.ndarray and len(val.shape)>1:
        varDim = ()
        for i,l in enumerate(val.shape):
            varDim_tmp = key + 'Sz' + str(i)
            varDim += (varDim_tmp,)
            ncid.createDimension(varDim_tmp,l)
        varType = get_type(val.flatten()[0])

    elif type(val) == np.ndarray or type(val) == list:
        varDim_tmp = key + 'Sz'
        varDim = (varDim_tmp,)
        ncid.createDimension(varDim_tmp,len(val))
        varType = get_type(val[0])
        # varType = val[0].dtype.char

    else:
        varDim = ('one',)
        varType = get_type(val)

    return varType, varDim


def compile_code():
    ## prepare code by compiling, if needed
    modes = '-g -Wall -ansi'

    libs = ' -lgsl -lgslcblas -lnetcdf -std=c++17' # or rather 11? (beta function requires 17)
    os.system('g++ %s -o ./theory/sharkfins ./theory/sharkfins.cpp %s' % (modes,libs))

    ## ---------------------------------------
    ## this is left just as an example of how to set paths (as I always forget...)
    #codepath = 'g++ -g -Wall -ansi  -o sharkfins -I/home/wollex/libs/gsl/include -I/home/wollex/libs/netcdf/include  sharkfins.cpp -L/home/wollex/libs/gsl/lib -lgsl -lgslcblas -L/home/wollex/libs/netcdf/lib -lnetcdf_c++ -lnetcdf -std=c++11'

    #libdir = ['gsl','netcdf']
    #pathlib = '/home/wollex/libs/'

    #buildL = ''
    #buildI = ''
    #for lib in libdir:
    #print lib
    #buildL += '-L' + pathlib + lib + '/lib '
    #buildI += '-I' + pathlib + lib + '/include '
    ## ---------------------------------------
