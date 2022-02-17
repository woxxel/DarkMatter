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

    def set_para(self,key,default,option,register=True):

        if key in option.keys():
            setattr(self,key,option[key])
        else:
            setattr(self,key,default)

        if register:
            self.paras.append(key)

    def prepare_sim_paras(self,steps):   # only works in conjunction with "set_simulation"
        # print(self.order)
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
    # defining model constants

    model = parameters('model')

    model.set_para('Npop',1,options)
    model.set_para('tau_M',0.010,options)        # membrane timeconstant in sec
    model.set_para('tau_A',0.005,options)        # GABAergic timeconstant in sec
    model.set_para('tau_N',0.2,options)          # NMDA timeconstant in sec
    model.set_para('J',-1.,options)              # synaptic strength

    model.set_para('kappa',1.,options)            # connectivity ratio (should have dim=Npop)
    model.set_para('drive',0,options)            # external drive (1=on, 0=off)

    if (model.drive > 0):
        model.set_para('tau_0',1.,options)         # external timeconstant
        model.set_para('J_0',1.,options)           # external synaptic weight
        model.set_para('K_0',1,options)            # external in-degree
    else:
        model.set_para('tau_0',0.,options)         # external timeconstant
        model.set_para('J_0',0.,options)           # external synaptic weight
        model.set_para('K_0',0,options)            # external in-degree

    ncid = Dataset(fileModel, "w");#, format="NETCDF4")
    ncid.createDimension('one',1)
    for key in model.paras:
        val = getattr(model,key)
        if type(val)==int:
            varType = np.dtype('int32').char
        else:
            varType = np.dtype('float64').char

        Var = ncid.createVariable(key,varType,('one',))
        Var[:] = val
    ncid.close()

    # model.print_parameter()

    return model

def set_simulation(fileSim,options,steps):
    # defining simulation parameters

    sim = parameters('simulation')

    sim.set_para('rateWnt',[0.,20.],options)
    sim.set_para('alpha_0',[0.,0.16],options)
    sim.set_para('tau_G',[0.005],options)
    sim.set_para('n',[0.],options)
    sim.set_para('eps',[0.],options)
    sim.set_para('eta',[0.],options)
    sim.set_para('I_alpha',[1],options)
    sim.set_para('I_beta',[1],options)

    sim.set_para('mode_calc',0,options)        # 0=exact, 1=approx
    sim.set_para('mode_stats',0,options)       # 0=general stats, 1=...

    # sim.set_para('nu0',0.,options)       #
    # sim.set_para('c',0.1,options)       #
    # sim.set_para('minZeta',-3.,options)       #
    # sim.set_para('maxZeta',3.,options)       #
    # sim.set_para('nZeta',11,options)       #

    order = ['rateWnt','alpha_0','tau_G','n','eps','eta','I_alpha','I_beta']
    sim.set_para('order',order,options)
    sim.prepare_sim_paras(steps)

    sv_str = 'mode=%d_steps=%d' % (sim.mode_stats,steps)
    sv_str += '_iter'
    for key in order:
        val = getattr(sim,key)
        if len(val) > 1:
            sv_str += '_%s=%g' % (key,val[-1])

    sv_str += '_const'
    for key in order:
        val = getattr(sim,key)
        if (len(val) == 1 or val[0]==val[-1]):
            sv_str += '_%s=%g' % (key,val[0])

    #write simulation parameters to netcdf file
    ncid = Dataset(fileSim, "w")
    ncid.createDimension('one',1)
    for key in sim.paras:
        val = getattr(sim,key)
        if type(val) == np.ndarray:
            varDim = key + 'Sz'
            ncid.createDimension(varDim,len(val))
            varType = val[0].dtype.char
        else:
            varDim = 'one'
            if type(val)==int:
                varType = 'i' #np.dtype('int32').char
            else:
                varType = 'f8' #np.dtype('float64').char
        if varType=='S':        ## somewhat complicated passing over readable strings to netcdf
            nChar = len(max(val,key=len))
            ncid.createDimension('charSz',nChar)
            Var = ncid.createVariable(key,'%s1'%(varType),(varDim,'charSz'))
            Var._Encoding = 'ascii'
            Var[:] = stringtochar(val)
        else:
            Var = ncid.createVariable(key,varType,(varDim,))
            Var[:] = val

    ncid.close()

    # sim.print_parameter()

    return sim, sv_str


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
