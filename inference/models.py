import numpy as np
import pymc3 as pm
import theano.tensor as tt

def construct_model_hierarchical(name,mP,mu_population=2.0,sigma_population=1.0,sigma_animal=1.0):

    """
        function to create a hierarchical model of the parameter 'name'.

        Input:
            name: string
                name of the parameter, used for distribution name
            mP: ModelParams Object
                including 'rates' DataFrame, arranged in a multicolumnar structure
                reflecting the different populations. Column names starting with '*'
                are used on the hierarchical level to infer mean value lower level
                distributions
            mu_population: float
                mean value of top-level prior
            sigma_population: float
                std of top-level prior
            sigma_animal:
                std of low-level prior

        Output:
            prior_animal: pm.distribution
                distribution containing a draw for each neuron in the dataset
                ! this should be masked, when populations contain different neuron
                numbers

        ToDo:
            * enable using different distributions, other than 'Normal'
    """

    # obtain information from data-column-header
    population_names = list(mP.rates.columns.names)
    population_shape = mP.get_data_shape()
    N_levels = len(population_names)
    # population_shapes = [
    #     len(np.unique(list(data.columns.get_level_values(name))))
    #     for name in population_names
    # ]

    # only columns with '*' are used for top-level inference
    is_hierarchical = [name.startswith('*') for name in population_names]

    hierarchical_shape = [
        sha if is_hierarchical[p] else 1
        for p,sha in enumerate(population_shapes)
    ]

    # create top-level distribution for prior mean-values
    prior_mean = mu_population + \
        pm.Normal(f'{name}_population',
            mu=0,
            sigma=sigma_population,
            shape=hierarchical_shape
        )
    tt.printing.Print('prior mean shape')(tt.shape(prior_mean))

    # create real distribution, from which values are drawn
    prior_animal = pm.Normal(f'{name}',
        mu=prior_mean,
        sigma=sigma_animal,
        shape=population_shapes
    )

    tt.printing.Print('prior animal')(tt.shape(prior_animal))
    prior_animal = prior_animal.reshape((-1,1))

    tt.printing.Print('prior animal')(tt.shape(prior_animal))

    # broadcast distribution to draw values for each neuron
    prior_animal = tt.tile(
        prior_animal,
        (mP.rates.shape[0], *[1]*N_levels)
    )

    return prior_animal

    # _num_animals, _num_layers, _num_clusters, nMax = shapes
    #
    # sigma_base = pm.HalfNormal(f'sigma_{name}',sigma=1.);
    # base = pm.HalfNormal(f'{name}_base',sigma=sigma_base);
    #
    # delta_population = pm.Normal(f'delta_{name}_layer',mu=0,sigma=sigma_population, shape=(1,_num_layers,_num_clusters))
    # #delta_cluster = pm.Normal(f'delta_{name}_cluster ',mu=0,sigma=sigma_cluster, shape=(1,1,_num_clusters))
    #
    # var = base + delta_population #+ delta_cluster
    #
    # animal = pm.Normal(f'{name}',mu=var,sigma=sigma_animal, shape=(_num_animals,_num_layers,_num_clusters))
    #
    # #pm.Deterministic(name,var)
    # animal = tt.tile(animal,(nMax,1,1,1),4)
    # tt.printing.Print(f'{name} shape')(tt.shape(animal))
    #
    # return animal#[mP.mask]


def construct_model_nu_max(name,shapes,nu_mean=20.,nu_sigma=5.):

    _num_animals, _num_layers, _num_clusters, nMax = shapes

    mu_nu_max = pm.Normal('mu_nu_max',mu=nu_mean,sigma=nu_sigma);
    sigma_nu_max = pm.Normal('sigma_nu_max',mu=5.,sigma=2.);

    # define nu_max as being able to differ between layers and clusters, but not animals
    nu_max = pm.Normal('nu_max',mu=mu_nu_max,sigma=sigma_nu_max,shape=(_num_animals,_num_layers,_num_clusters));
    nu_max = tt.tile(nu_max,(nMax,1,1,1),4)

    tt.printing.Print('nu')(tt.shape(nu_max))

    return nu_max#[mP.mask]
