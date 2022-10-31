import pymc3 as pm
import theano.tensor as tt

def construct_model_hierarchical(name,mP,sigma_animal=1.0,sigma_layer=1.0,sigma_cluster=1.0):

    ## build gamma as a base value + modifications per layer and per animal
    # gamma = gamma_base + delta_gamma_layer (N(0,sigma_layer)) + delta_gamma_animal (N(0,sigma_animal))

    sigma_base = pm.HalfNormal(f'sigma_{name}',sigma=1.);
    base = pm.HalfNormal(f'{name}_base',sigma=sigma_base);

    delta_layer = pm.Normal(f'delta_{name}_layer',mu=0,sigma=sigma_layer, shape=(1,mP._num_layers,1))
    #delta_cluster = pm.Normal(f'delta_{name}_cluster',mu=0,sigma=sigma_cluster, shape=(1,1,mP._num_clusters))

    var = base + delta_layer #+ delta_cluster

    animal = pm.Normal(f'{name}',mu=var,sigma=sigma_animal, shape=(mP._num_animals,mP._num_layers,mP._num_clusters))

    #pm.Deterministic(name,var)
    animal = tt.tile(animal,(mP.nMax,1,1,1),4)
    tt.printing.Print(f'{name} shape')(tt.shape(animal))

    return animal[mP.mask]


def construct_model_nu_max(name,mP,nu_mean=20.,nu_sigma=5.):

    mu_nu_max = pm.Normal('mu_nu_max',mu=nu_mean,sigma=nu_sigma);
    sigma_nu_max = pm.Normal('sigma_nu_max',mu=5.,sigma=2.);

    # define nu_max as being able to differ between layers and clusters, but not animals
    nu_max = pm.Normal('nu_max',mu=mu_nu_max,sigma=sigma_nu_max,shape=(mP._num_animals,mP._num_layers,mP._num_clusters));
    nu_max = tt.tile(nu_max,(mP.nMax,1,1,1),4)

    tt.printing.Print('nu')(tt.shape(nu_max))

    return nu_max[mP.mask]
