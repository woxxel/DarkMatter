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

    def set_para(self,key,default,option,register=True,sz=None):

        val = option[key] if key in option.keys() else default
        if sz:
            if not type(val)==list:
                val *= np.ones(sz)
            else:
                assert len(val)==sz, f'Please provide values for {key} of dimension 1 or {sz}'

        setattr(self,key,val)
        # if key in option.keys():
        # else:
        #     setattr(self,key,default)

        if register:
            self.paras.append(key)

    def prepare_sim_paras(self,steps):   # only works in conjunction with "set_simulation"
        for key in self.order:
            val = getattr(self,key)
            assert type(val) == list, 'Please specify all iteration parameters as lists!'
            if len(val) == 2:
                val = np.linspace(val[0],val[-1],steps+1)[1:]#.astype('float')
            elif len(val) > 2:
                val = np.array(val)
            else:
                val = np.array(val).astype('float')
            setattr(self,key,val)

        nChar = len(max(self.order,key=len))
        self.order = np.array(self.order,dtype='S'+str(nChar))

    def print_parameter(self):

        for key in self.paras:
            print(key, getattr(self,key))

def set_model(fileModel,options):

    # defining model constants (all of them!)
    model = parameters('model')

    ## global parameters
    model.set_para('Npop',2,options)
    Npop = getattr(model,'Npop')
    assert not (Npop % 2), f'Number of populations must be an even number, Npop={Npop} given'

    model.set_para('tau_M',0.010,options)        # membrane timeconstant in sec

    ## set arbitrary number of synaptic timeconstants
    model.set_para('tau_order',[0,0,1],options)
    model.set_para('tau_I',[0.005,0.2,0.03],options)
    model.set_para('tau_n',[1.,0.,1.],options)
    model.set_para('tau_norm',[1.,1.,1.],options)

    assert not len(set(range(Npop)) - set(getattr(model,'tau_order'))), 'Please make sure there are synaptic timeconstant parameters defined for each population'

    for tau_para in ['tau_I','tau_n','tau_norm']:
        assert len(getattr(model,tau_para)) == len(getattr(model,'tau_order')), 'Make sure to provide the same number of parameters for each variable of synaptic timeconstants'
        assert len(getattr(model,tau_para)) >= Npop, 'Number of defined synaptic timeconstants needs to be at least as large as the number of populations'
    # model.set_para('tau_A',0.005,options)        # GABAergic timeconstant in sec
    # model.set_para('tau_N',0.2,options)          # NMDA timeconstant in sec

    # print(Npop)
    ### parameters for pairs of populations
    model.set_para('eps',0.,options,sz=Npop//2)
    model.set_para('eta',0.,options,sz=Npop//2)

    ### parameters for each population
    model.set_para('J',-1.,options,sz=Npop)              # synaptic strength
    model.set_para('I_ext',1,options,sz=Npop)           # external, constant drive
    model.set_para('kappa',1.,options,sz=Npop)            # connectivity ratio (should have dim=Npop)
    model.set_para('drive',0,options,sz=Npop)            # external drive (1=on, 0=off)
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
        # if type(val)==int:
        #     varType = np.dtype('int32').char
        # else:
        #     varType = np.dtype('float64').char

        Var = ncid.createVariable(key,varType,varDim)
        Var[:] = val
    ncid.close()

    model.print_parameter()

    return model

def set_simulation(fileSim,options,steps):

    # defining simulation parameters (only the ones, which are being iterated)
    sim = parameters('simulation')

    for key in options['simulation']:
        if not (key in ['sim_p','sim_tau']):
            sim.set_para(key,[],options['simulation'])

    sim.set_para('order',sim.paras,{},register=False)

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
    for key in sim.paras:

        val = getattr(sim,key)
        varType,varDim = prepare_ncid_var(ncid,sim,key,val)

        if varType=='S':        ## somewhat complicated passing over readable strings to netcdf
            nChar = len(max(val,key=len))
            ncid.createDimension('charSz',nChar)
            Var = ncid.createVariable(key,'%s1'%(varType),varDim+'charSz')
            Var._Encoding = 'ascii'
            Var[:] = stringtochar(val)
        else:
            Var = ncid.createVariable(key,varType,varDim)
            Var[:] = val

    ncid.close()

    # sim.print_parameter()
    # print(sim.paras)

    return sim, sv_str

def prepare_ncid_var(ncid,sim,key,val):

    if type(val) == np.ndarray or type(val) == list:
        varDim = key + 'Sz'
        ncid.createDimension(varDim,len(val))

        if type(val[0])==int:
            varType = np.dtype('int32').char #'i' #
        else:
            varType = np.dtype('float64').char
        # varType = val[0].dtype.char
    else:
        varDim = 'one'
        if type(val)==int:
            varType = np.dtype('int32').char #'i' #
        else:
            varType = np.dtype('float64').char #'f8' #

    return varType, (varDim,)


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
