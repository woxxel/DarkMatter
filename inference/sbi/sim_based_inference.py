from scipy.stats import ks_2samp

import torch
from sbi import utils as utils
from sbi import analysis as analysis
from sbi.inference.base import infer

from empirical.read_data import *
from darkMatter.darkMatter import *

# create class with simulator method for sbi
class sbi_darkmatter:

    data = None
    animal = 2
    dataPath = '../data/BuscheLab/spiking_data_for_modeling.mat'
    N = None
    T = None

    plot = False

    def __init__(self,plot=False):

        self.plot = plot


    def register_data(self):

        # read in data from file and prepare cumulative distribution
        self.data = read_data(self.dataPath)

    def create_cdf(self,animal=1):

        # this is not the animal number, but the data-set-number. fix this!
        data = self.data[animal]['rate'][:,0]

        # create cdf
        self.N = len(data)
        self.cdf_x = np.sort(data)
        self.cdf_y = np.arange(self.N) / float(self.N)

    def create_prior(self):
        # uniform priors on all?
        # maybe rather Gaussian/student-t?

        """
            alpha_0: uniform
            tau_I (3x): student-t around empirical data
            tau_n: gauss around empirical?
            rate: student-t around measured stuff
        """

        pass


    def simulator(self,param_set):
        """
            param_set: list(float)
                contains parameters in the order [rate,alpha_0,tau_I,tau_n]
                but actually: should/could have gamma/delta as input and subsequently find parameters
        """

        # should call c++ with as little overhead as possible (what can be avoided? what not?)
        # print('params: ',param_set)
        rate, alpha_0, tau_G = param_set

        L = 1
        J_l = np.ones((L,L))
        np.fill_diagonal(J_l,0)

        options = {
            # count of layers, populations, PSPs
            'L': L,
            'P': 2,
            'S': [1,2],     # contains number of synapses for each population

            # layer level parameters
            'eps': 0,
            'eta': 0.9,
            'J0_l': J_l,
            'kappa': 1.,

            # population level parameters
            'I_ext': 1,
            'rateWnt': rate,
            'alpha_0': alpha_0,
            'tau_M': 0.01,
            'J0': 1.,

            # psp level parameters
            'tau_I': [tau_G,0.005,0.2],
            'tau_n': 0.,
            'tau_norm': 1.,

            # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
            'mode': 1,      # mode 0: phase plots, mode 1: data simulation
            'mode_stats': 0,
            'mode_calc': 0,
            'computation': {
                'N': self.N,    # number of neurons to draw from
                'T': 600.,   # time of measurement

                'draw_from_theory': 100,
                'draw_finite_time': 1,
                # 'seed': 10,
            }
        }

        res = darkMatter(steps=1,options=options,mode=1,rerun=True)

        pvals = np.zeros(options['computation']['draw_from_theory'])
        for i in range(res['rates_T'].shape[-2]):
            dat = ks_2samp(res['rates_T'][:,i,0],self.cdf_x)
            pvals[i] = dat.pvalue

        if self.plot:
            fig,ax = plt.subplots(2,2,figsize=(8,4))
            ax[0][1].hist(self.cdf_x,bins=np.linspace(0,max(self.cdf_x),51),alpha=0.5,density=True)

            ax[0][0].plot(self.cdf_x,self.cdf_y)

            worst = [None,0]
            best = [None,1]
            for i in range(res['rates_T'].shape[-2]):
                N = res['rates_T'].shape[0]
                x = np.sort(res['rates_T'][:,i,0])
                y = np.arange(N) / float(N)
                ax[0][0].plot(x,y,lw=0.5,color='gray')

            best = np.argmin(pvals)
            worst = np.argmax(pvals)

            ax[0][1].hist(res['rates_T'][:,worst,:].flatten(),bins=np.linspace(0,max(self.cdf_x),51),alpha=0.5,density=True,color='red')
            ax[0][1].hist(res['rates_T'][:,best,:].flatten(),bins=np.linspace(0,max(self.cdf_x),51),alpha=0.5,density=True,color='green')

            ax[1][0].hist(pvals,np.logspace(-4,0,26))
            ax[1][0].set_xscale('log')

            plt.show(block=False)

        # print(f'pvalue {pvals.mean()} +/- {pvals.std()}')
        return_val = np.zeros((1,1))
        return_val[0,0] = pvals.mean()
        # return_val = np.zeros((1,2))
        # return_val[0,:] = [pvals.mean(),pvals.std()]
        return return_val

        # return res

SBI = sbi_darkmatter()
SBI.register_data()
SBI.create_cdf(animal=2)

prior = utils.BoxUniform(low=np.zeros(3), high=np.array([10,0.2,0.2]))
posterior = infer(SBI.simulator,prior, method='SNPE',num_simulations=1000)

observation = torch.ones(1)

samples = posterior.sample((100000,), x=observation)
log_probability = posterior.log_prob(samples, x=observation)
_ = analysis.pairplot(samples, limits=[[0,10],[0,0.2],[0,0.1]], figsize=(12,12))
plt.show(block=False)
