import numpy as np

def set_options(L=1,nI=1,nE=1):
    '''
        this function sets the parameters for the simulation
        it contains default values for all

        TODO:
            * add proper description for all parameters
            * should translate to class, in which primary & secondary can be set 
                by methods?
            * add kwargs to set non-default values
    '''
    
    ## setting network parameters
    kappa = create_population_values(1/nI if nI else 0,nI,1/nE if nE else 0,nE)
    J0 = create_population_values(-1.,nI,1.,nE)
    S = create_population_values(1,nI,2,nE)
    tau_I = create_population_values(0.01,nI,[0.005,0.2],nE)
    
    J_l = np.ones((L,L))
    np.fill_diagonal(J_l,0)

    options = {
        # count of layers, populations, PSPs
        'L': L,
        'P': len(S),
        'S': S,     # contains number of synapses for each population

        # layer level parameters
        'eps': 1./np.sqrt(2),
        'eta': 0.9,
        'J0_l': J_l,

        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'kappa': kappa,
        'alpha_0': 0.01,
        'Psi_0': 0.,
        'tau_M': 0.01,
        'J0': J0,

        # psp level parameters
        'tau_I': tau_I,
        'tau_n': 0.3,
        'tau_norm': 1.,

        'mode': 0,
        'mode_stats': 0,    ## stats:
        ###    0: sharkfins
        ###    1: rate_dependence stats
        ###    2: borders phase space
        ###    2: compare exact vs. approx (single)
        ###    3: KL-phase-space (costly!)
        
        'mode_calc': 0,
        'mode_selfcon': 0,  
        ### 0 = mean rates for all populations provided, 
        ### 1 = mean rate across populations provided, assuming same alpha

        'simulation': {
            ## this dictionary contains the extents of the parameters space to be
            ## sampled and has to be set before running the calculations

            # 'Psi_0': [-0.15,0.15],     # first parameter set is "primary" one
            # 'alpha_0': [0.,0.02],      # the other one is the "secondary"
            
            # for each iteration parameter, specify (layer,population,psp)-tuple
            # specify -1 if should be applied to all candidates
            # 'sim_prim': [0,1,0],       # when population parameters are iterated, specify population number(s) (negative = all)
            # 'sim_sec': [0,-1,0],
        }
    }
    return options


def create_population_values(I_val,nI,E_val,nE):

    I_val = I_val if type(I_val)==list else [I_val]
    E_val = E_val if type(E_val)==list else [E_val]

    val = I_val*nI
    val.extend(E_val*nE)
    return val


def get_J(eta,eps,tau_M,p_pre,p_post):

    if p_post=='E' and p_pre=='E':
        return eta*eps*tau_M
    elif p_post=='E' and p_pre=='I':
        return np.sqrt(1-(eta*eps)**2)*tau_M
    elif p_post=='I' and p_pre=='E':
        return eps*tau_M
    elif p_post=='I' and p_pre=='I':
        return np.sqrt(1-eps**2)*tau_M


def get_var_V(J_from_E,J_from_I,rateWnt,tau_A,tau_N,tau_G,tau_M,n):
    ## variance of the membrane potential
    ## from excitatory AMPA synapses

    var = {}
    var['V_A'] = J_from_E**2 * rateWnt / (tau_A + tau_M) * ( (1-n)**2/2 + (1-n)*n*tau_A / (tau_A + tau_N) )
    ## from excitatory NMDA synapses
    var['V_N'] = J_from_E**2 * rateWnt / (tau_N + tau_M) * (n**2/2 + (1-n)*n*tau_N / (tau_A + tau_N) )
    ## from inhibitory GABA synapses
    var['V_G'] = J_from_I**2 * rateWnt * 0.5 / (tau_G + tau_M)

    N = 1
    for key in var.keys():
        N = max(N,var[key].shape[0] if isinstance(var[key],np.ndarray) else 1)
    
    for key in var.keys():
        if not isinstance(var[key],np.ndarray):
            var[key] = np.full(N,var[key])
        # print(var)
    
    return var['V_A'],var['V_N'],var['V_G']
