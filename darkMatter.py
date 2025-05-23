import numpy as np
import os
from pathlib import Path
# from utils.parameters import parameters
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

    TODO:
        * expand explanation for "options" parameter
        * add explanation for "mode" parameter (needed?)
        * add overall comment (written by, at, ...)
"""

def darkMatter(steps=100,mode=0,options={},path=None,logging=3,cleanup=True,rerun=False,compile=False,suffix=''):

    if suffix != '':
        suffix = f'_{suffix}' if not suffix.startswith('_') else suffix

    print('DM suffix:',suffix)
    # set model parameters
    path_script = os.path.dirname(__file__)
    path = path_script if (path is None) else path

    # print('script:',path_script)
    # print('path:',path)

    fileModel = Path(path,'data',f'modPara{suffix}.nc')

    model,sv_str_const = set_model(fileModel,options)
    if mode==0:
        # print('hello')
        filePara = Path(path,'data',f'simPara{suffix}.nc')
        sim, sv_str = set_simulation(filePara,options,steps)
    elif mode==1:
        filePara = Path(path,'data',f'comPara{suffix}.nc')
        com, sv_str = set_computation(filePara,options['computation'])

    # fileResults = './.data/results.nc'
    fileResults = Path(path,'data',f'results_{sv_str}_{sv_str_const}{suffix}.nc')
    if not fileResults.is_file() or rerun:
        if mode==0:
            program = Path(path_script,'theory','sharkfins')
        elif mode==1:
            program = Path(path_script,'theory','generate_measures')

        if (compile):
            compile_code(path_script,program)

        run_str = f'{program} {fileModel} {filePara} {fileResults} {logging}'
        os.system(run_str)

    # obtain data from simulation
    results = read_from_netcdf(fileResults)

    if cleanup:
        for f in [fileModel,filePara,fileResults]:
            os.remove(f)
    else:
        results['filePaths'] = [fileModel,filePara,fileResults]

    return results


def read_from_netcdf(path):

    ncid = Dataset(path, "r")

    data = {}
    for key in ncid.variables:
        data[key] = ncid.variables[key][:]
    ncid.close()

    return data


def write_to_netcdf(class_obj, path):
    # write simulation parameters to netcdf file
    ncid = Dataset(path, "w")
    try:
        ncid.createDimension("one", 1)
        # sim.paras.extend(['mode_calc','mode_stats','mode_selfcon','sim_prim','sim_sec'])#,'order'
        for key in class_obj.paras:

            val = getattr(class_obj, key)
            varType, varDim = prepare_ncid_var(ncid, key, val)
            if varType == "S":
                nChar = len(max(val, key=len))
                ncid.createDimension("charSz", nChar)
                Var = ncid.createVariable(key, "S1", varDim + ("charSz",))
                Var[:] = np.array([list(s.ljust(nChar)) for s in val], "S1")
            else:
                Var = ncid.createVariable(key, varType, varDim)
                Var[:] = val
        ncid.close()
    except:
        ncid.close()
        assert f"Error in writing parameters to {path}"


def set_model(fileModel,options):

    # defining model constants (all of them!)
    model = parameters('model')

    ## global parameters
    model.set_para("L", 1, options, np.int16, sz=1)
    L = np.sum(getattr(model,'L'))

    # print(f"{model.L=}")
    # print(f"{model.L.size=}")
    # print(len(model.L))

    model.set_para("P", 2, options, np.int16, sz=L)
    P = getattr(model,'P')       # populations per layer
    nP = np.sum(P)    # total number of populations

    model.set_para("S", [1, 2], options, np.int16, sz=nP)
    S = getattr(model,'S')
    nS = np.sum(S)

    ### parameters for layers
    model.set_para('eps',0.,options,sz=L)
    model.set_para('eta',0.,options,sz=L)
    model.set_para('J0_l',1.,options,sz=L)              # synaptic strength of inter-layer interactions

    ### parameters for each population
    model.set_para(
        "I_ext", -1.0, options, dtype=float, sz=nP
    )  # external, constant drive
    model.set_para('rateWnt',1.,options,sz=nP)           # external, constant drive
    model.set_para('kappa',1.,options,sz=nP,register=False)            # connectivity ratio (should have dim=Npop)
    model.set_para('alpha_0',0.,options,sz=nP)           # external, constant drive
    model.set_para('Psi_0',0.,options,sz=nP,register=False)           # external, constant drive
    model.set_para('tau_M',0.01,options,sz=nP,register=False)        # membrane timeconstant in sec
    model.set_para('tau_n',0.,options,sz=nP)
    model.set_para('J0',-1.,options,sz=nP,register=False)             # base synaptic strength

    model.set_para(
        "drive", 0, options, np.int16, sz=L, register=False
    )  # external drive (1=on, 0=off)

    ### parameters for each synapse set
    model.set_para('tau_I',[0.005,0.2,0.03],options,sz=nS)
    model.set_para('tau_norm',1.,options,sz=nS,register=False)

    for tau_para in ['tau_I','tau_norm']:
        assert len(getattr(model,tau_para)) == nS, 'Make sure to provide the same number of parameters for each variable of synaptic timeconstants'
        assert len(getattr(model,tau_para)) >= nP, 'Number of defined synaptic timeconstants needs to be at least as large as the number of populations'

    ## write parameters for external spiking input (not currently implemented)
    # for p in range(Npop):
    #     if (model.drive[p] > 0):
    #         model.set_para('tau_0',1.,options,sz=Npop)         # external timeconstant
    #         model.set_para('J_0',1.,options,sz=Npop)           # external synaptic weight
    #         model.set_para('K_0',1,options,sz=Npop)            # external in-degree
    #     else:
    #         model.set_para('tau_0',0.,options,sz=Npop)         # external timeconstant
    #         model.set_para('J_0',0.,options,sz=Npop)           # external synaptic weight
    #         model.set_para('K_0',0,options,sz=Npop)            # external in-degree

    model.set_para('I_alpha',1.,options,register=False)
    model.set_para('I_beta',1.,options,register=False)

    write_to_netcdf(model, fileModel)

    sv_str = '_const'
    for i,key in enumerate(model.paras):
        if model.register[i]:
            val = getattr(model,key)
            # print(key, val)
            # if np.isscalar(val) and isinstance(val,int):#or
            if val.size == 1 and isinstance(val, int):  # or
                # print("scalar int")
                sv_str += f"_{key}={val.flat[0]:d}"
            # elif np.isscalar(val):
            elif val.size == 1:
                # print("scalar non-int")
                sv_str += f"_{key}={val.flat[0]:.3g}"
            elif (val.size == 1 or val[0] == val[-1]) and isinstance(val[0], float):
                # print("same vals, float")
                sv_str += f'_{key}={val[0]:.3g}'
            elif (
                val.size == 1 or val[0] == val[-1]
            ):  # and (isinstance(val[0],int) or val[0]%1==0):
                # print("same vals, else")
                sv_str += f'_{key}={val[0]:d}'

    # model.print_parameter()

    return model, sv_str


def set_simulation(fileSim,options,steps):

    # defining simulation parameters (only the ones, which are being iterated)
    sim = parameters('simulation')

    for key in options['simulation']:
        if key in ['sim_prim','sim_sec']:
            sim.set_para(key, [], options["simulation"], np.int16, register=False)
        else:
            sim.set_para(key, [], options["simulation"], np.int16)

    # sim.set_para('order',sim.paras,{},register=False)

    # sim.set_para('I_alpha',[1],options)
    # sim.set_para('I_beta',[1],options)

    sim.prepare_sim_paras(steps)

    sim.set_para(
        "mode", 0, options, np.int16, register=False
    )  # 0=phase plots, 1=data simulation
    sim.set_para("mode_calc", 0, options, np.int16, register=False)  # 0=exact, 1=approx
    sim.set_para(
        "mode_stats", 0, options, np.int16, register=False
    )  # 0=general stats, 1=...
    sim.set_para(
        "mode_selfcon", 0, options, np.int16, register=False
    )  # 0=obtain 2nd moments, 1=obtain 2nd moments and all but first rate

    # order = ['rateWnt','alpha_0','tau_G','n','eps','eta','I_alpha','I_beta']
    order = [o for o in list(options["simulation"]) if not o.startswith("sim")]
    sim.set_para("order", order, options, register=False)

    write_to_netcdf(sim, fileSim)

    sv_str = (
        f"mode={sim.mode}_stats={sim.mode_stats}_approx={sim.mode_calc}_steps={steps}"
    )
    sv_str += '_iter'
    for i,key in enumerate(sim.paras):
        if sim.register[i]:
            val = getattr(sim,key)
            if len(val) > 1:
                # sv_str += '_%s=%g' %^(key,val[-1])
                sv_str += f'_{key}'

    # sim.print_parameter()

    return sim, sv_str


def set_computation(fileComputation,options):

    # defining model constants (all of them!)
    com = parameters('computation')

    ### parameters for layers
    com.set_para("N", 10, options, np.int16)
    com.set_para('T',100,options)
    com.set_para("draw_from_theory", 100, options, np.int16, register=False)
    com.set_para("draw_finite_time", 100, options, np.int16, register=False)
    # com.set_para('seed_theory',np.random.randint(0,10000,getattr(com,'draw_finite_time')),options)

    write_to_netcdf(com, fileComputation)

    sv_str = 'computation'
    for i,key in enumerate(com.paras):
        if com.register[i]:
            val = getattr(com,key)
            if (type(val)!=list and type(val)!=np.ndarray):
                sv_str += f'_{key}={val:.0f}'

    com.set_para("seed", np.random.randint(0, 2**16), options, np.int16)

    # com.print_parameter()

    return com, sv_str


class parameters:

    # class to support writing datafiles
    def __init__(self,name):
        # print('Building %s class'%name)
        self.paras = []
        self.register = []

    def set_para(
        self, key, default, option=None, dtype=np.float32, register=True, sz=None
    ):

        if option:
            val = option[key] if key in option.keys() else default
        else:
            val = default

        # print(key, val)
        # print(type(val), np.array(val).dtype)
        # if not (type(val) == list and type(val[0]) == str):
        #     val = np.array(val).astype(dtype)

        if sz:
            if not (type(val)==list or type(val)==np.ndarray):
                val = np.full(sz,val)
                # val *= np.ones(sz)
            else:
                if np.array(val).size == 1:
                    val = np.full(sz, val)
                else:
                    f = sz//len(val)
                    val = np.array(list(val)*f)
                assert (
                    np.array(val).size == sz
                ), f"Please provide values for {key} of dimension 1 or {sz}"

        setattr(self,key,val)
        # if key in option.keys():
        # else:
        #     setattr(self,key,default)

        self.paras.append(key)
        self.register.append(register)

    def prepare_sim_paras(self,steps):   # only works in conjunction with "set_simulation"
        for key in self.paras:
            val = getattr(self,key)
            print(key, val, type(val))
            assert (
                type(val) == list or type(val) == np.ndarray
            ), f"{key} not provided as a list ({val}). Please specify all iteration parameters as lists!"
            if len(val) == 2:
                val = np.linspace(val[0],val[-1],steps+1)[1:]#.astype('float')
            elif len(val) > 2:
                val = np.array(val)
            else:
                val = np.array(val).astype('float')
            setattr(self,key,val)

        nChar = len(max(self.paras,key=len))
        # print("order:", self.paras)
        self.order = np.array(self.paras,dtype='S'+str(nChar))

    def print_parameter(self):

        for key in self.paras:
            print(key, getattr(self,key))


def get_type(val):
    if type(val)==np.int_ or type(val)==int:
        varType = np.dtype('int32').char #'i' #
    elif type(val) == np.bytes_ or type(val) == bytes or type(val) == str:
        varType='S'
    else:
        varType = np.dtype('float64').char
    return varType


def prepare_ncid_var(ncid,key,val):

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


def compile_code(path_script,program):
    ## prepare code by compiling, if needed
    
    path = Path(path_script,'theory','src','3rdParty')
    
    modes = '-g -Wall -ansi'
    manual_libs = f"-I{path / 'include'} -L{path / 'lib'}"
    module_libs = ' -lgsl -lgslcblas -lnetcdf -std=c++17' # or rather 11? (beta function requires 17)

    cmd = f"g++ {modes} -o {program} {program}.cpp {manual_libs} {module_libs}"

    print(cmd)
    os.system(cmd)

    ## ---------------------------------------
    ## this is left just as an example of how to set paths (as I always forget...)

    #libdir = ['gsl','netcdf']
    #pathlib = '/home/wollex/libs'

    #buildL = ''
    #buildI = ''
    #for lib in libdir:
    #   print lib
    #   buildL += os.path.join('-L',pathlib,lib,'lib ')
    #   buildI += os.path.join('-I',pathlib,lib,'include ')

    #codepath = 'g++ -g -Wall -ansi  -o sharkfins ' + \
    #  buildI + ' sharkfins.cpp ' + \
    #  buildL + ' -lgsl -lgslcblas -lnetcdf_c++ -lnetcdf -std=c++11'

    ## ---------------------------------------
