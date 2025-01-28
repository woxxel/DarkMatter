#include <vector>
// #include <random>
#include "gsl/gsl_rng.h"
#include <iostream>

using namespace std;

class Results;
class Simulation;
class Model;


struct Model_Results
{
    vector< vector<int> > trans_inc, trans_imp;
    // vector< vector<double> > trans_inc, trans_imp;

    void initiate(unsigned sim_steps_1, unsigned sim_steps_2, unsigned sim_mode);
};

struct Population_Results
{
    vector<vector<double> > rate, q, alpha_raw, alpha, sigma_V, gamma, delta, rate_max, chi, I_balance, regions, implausible, entropy, KL_entropy;


    vector<vector<double> > infoContent;

    // vector<unsigned> axesDim;

    // vector<vector<double> > KS;
    // vector<double> KL_entropy;

    // vector<vector<vector<double> > > rate_T;
    // vector<vector<vector<int> > > N_AP;
    // vector<vector<double> > rate_inf;
    // vector<double> p_est_inf, p_est_T;
    // vector<vector<vector<double> > > p_bayes_est;
    // vector<double> p_bayes_est_measures;
    // vector<double> cdf_theory;
    // vector<double> p_k;
    // vector<int> p_hist;


    vector< vector<int> > trans_DM, trans_np;
    // vector< vector<double> > trans_DM_approx, trans_np_approx;

    // double d_nu, max_prob, factor;
    // int steps;

    int axes[7] = {0,0,0,0,0,0,0};

    void initiate(unsigned sim_steps_1, unsigned sim_steps_2, unsigned sim_mode);
};

struct Population_Simulation
{
    bool in_DM, in_np;
    int trans_DM, trans_np;
    bool trans_DM_found = false, trans_np_found = false;

    double rate, q;
    double sigma_V, alpha, alpha_raw;
    double I_balance, gamma, delta, rate_max, chi;
    double regions, implausible;
    double infoContent, KL, entropy;
    double max_prob;

    double distribution_exact(double nu);

    double draw_rate(gsl_rng *rng);
    double draw_sample(gsl_rng *rng, double rate, double T);

};

struct Model_Simulation
{
    bool in_inc, in_imp;
    int trans_inc, trans_imp;
    bool trans_inc_found = false, trans_imp_found = false;
    size_t nTrans = 0;
};

struct PSP
{
    /*
        Definition of synaptic timeconstants
            tau_I: values in secs
            tau_norm: normalization
                <0: abs(tau_norm) = value of peak height
                >0: tau_norm = total current transported (integral)
                =0: invalid
            tau_n: mix-ratio - has to be normalized s.t. sum(tau_n) = 1
    */
    unsigned s; // PSP ID
    double tau_I, tau_norm, tau_n;
    void print_PSP();
};

struct Population
{
    // Defines parameters for each population of the network.
    // Each excitatory population is accompanied by an inhibitory pair, whos interactions are specified by eta & eps

    unsigned p; // population ID
    vector<PSP> psp;
    unsigned nPSP;

    Population_Results results;
    Population_Simulation simulation;

    /*
        External drive to populations
            I_ext (vector<bool>)
                true: introduces external drive to keep population firing rate at rate specified by rateWnt
                false: no external drive - population is only activated by excitatory populations
            rateWnt (vector<double>)
                if I_ext[p], rateWnt defines the rate at which the population should be active - otherwise neglected

            drive (vector<double>)
                specifies (along with 2-3 other variables) spiking, external drive, introducing further heterogeneity into the network
    */
    // parameters:
	//  	drive:
	//  	0 : only excitatory one driven
  //  	1 : all populations driven to have same firing rates
	//  	2 : recurrent, inhibitory population is driven by afferent spikes from excitatory population

    int I_ext;
    double rateWnt;
    double tau_n;
    double J0; // synaptic strength base value
    double kappa;

    // int drive;
    // double K_0; // average incoming number of synapses as multiple of recurrent connections number K
    // double tau_0;

    /*
        further parameters
            alpha_0 (double)
                heterogeneity of this population
    */
    double alpha_0, Psi_0;

    double tau_M;
    vector<double> J;

    void print_population();

};

struct Layer
{
    // Defines parameters for each layer of the network.
    // Each layer is composed of an excitatory and an inhibitory population, whos interactions are specified by eta & eps
    vector<Population> population;
    unsigned nPop;
    unsigned l; // layer ID


    /*
        Interaction between inhibitory and excitatory population
            eps (vector<double>)
                strength of excitatory-inhibitory feedback loop (0 < eps < 1/sqrt(2))
            eta (vector<double>)
                inter-population excitatory coupling (0 < eta < 1)

            J_l (vector<double>)
                specifies the strength of incoming synapses from layer l (only to excitatory for now)
            kappa (vector<double>)
                specifies the relative number of inhibitory neurons to excitatory ones
    */
    double eta, eps;
    vector<double> J0_l, J_l;

    void print_layer();
};

struct parameters
{
    unsigned Npop;

    // vector<population> pop; // holds all populations of the network

    // vector<int> tau_order;

	double tau_G, tau_A, tau_N, tau_M, J;
	double eta, eps, n;
    double kappa;
    vector<double> rate, q, alpha_0, Psi_0;

    //external drive
    // int drive;
    // double J0; // synaptic strength with respect to recurrent weights J
    // double K_0; // average incoming number of synapses as multiple of recurrent connections number K
    // double tau_0;

	vector<double> J_E, J_I;
//         vector<double> kappa;

	vector<double> alpha_raw, alpha, sigma_V;
	vector<double> rate_max;
	vector<double> gamma, delta, chi, I_balance;
    vector<double> regions;
    vector<double> entropy, KL;
};


struct info_paras
{
    double nu0, c;
};

class Model
{
	public:
        unsigned L;
        vector <Layer> layer;

        unsigned nPop;
        double network_rate; // needed, e.g. in the case of split subpopulations

        double I_alpha, I_beta;

        Model_Simulation simulation;
        Model_Results results;

        parameters paras;

        vector<double> infoContent;

        // void add_PSP(int p, double tau_I, double tau_norm, double tau_n);
        void set_parameters();
        void set_rates();
        void set_weights();
        void set_mixture();
        void solve_selfcon(int mode_calc);
        void solve_selfcon_split(int mode_calc);
        void solve_selfcon_from_currents(int mode_calc);
        void write_results();

        bool q_border1(Population_Simulation *popSimP);
        bool q_border2(Population_Simulation *popSimP);

        bool is_darkMatter(Population_Simulation *popSimP);
        bool is_noPeak(Population_Simulation *popSimP);
        bool is_inconsistent(Model_Simulation *mSimP);
        bool is_implausible(Model_Simulation *mSimP);

        void integrate_information();

        void find_transitions(Simulation *simP);
        void store_update();


// 		model();
        vector< vector<double> > calc_alpha(vector< vector<double> > q);
        // void resize();

        void get_max_prob();
        void get_sigma_V();
        
        double nu_peak_log_full(Population_Simulation *simP);

    private:

        void get_alpha();
        void get_delta();
        void get_gamma();
        void get_chi();

};

struct Simulation_Variable
{
    double *valP, *paraP_approx;
    vector<double*> paraP;

    unsigned paraSz;

    string name;
    unsigned iter, steps;

    // void (model::*onUpdate)();    // function which is called every time the variable updates

    void initialize(double *modVarP, double *varP, unsigned varSz, string varName);
    bool iterate();
    void print_status();

};


class Simulation
{
    public:

        info_paras infoParas;

        // bool in_imp, in_inc;
        bool trans_imp_found, trans_inc_found, trans_imp_found_approx, trans_inc_found_approx;
        double trans_inc, trans_imp;

        vector<double> eps, eta, alpha_0, Psi_0, rateWnt, I_alpha, I_beta, tau_I, tau_n;
        vector<vector<int> > sim_pointer;

        unsigned nVar;
        vector<string> order; // contains the strings of parameters in the order that the iteration should walk through

        // unsigned n_iter, alpha_0_iter, tau_G_iter, rateWnt_iter, eps_iter, eta_iter;
        int mode_calc, mode_stats, mode_selfcon;
        size_t tau_nSz, alpha_0Sz, Psi_0Sz, tau_ISz, rateWntSz, epsSz, etaSz, I_alphaSz, I_betaSz, orderSz, charSz, steps, sim_primSz, sim_secSz;

        // vector<bool> trans_DM_found, trans_np_found, trans_DM_found_approx, trans_np_found_approx;

        double y_val();

        vector<Simulation_Variable> vars; // initialize, such that variables are registered in order
        bool iter_status[6];
        // bool iterating = true;
        void set_vars(Model *modP, unsigned idx, unsigned l, unsigned p, unsigned s);
        void initialize(Model *modP);
        void reset_iteration();
        bool run_iteration(Model *modP, Model *modP_approx);
        void print_simStatus();

        void initiate_y_axis(Model *modP);
        void store_results(Model *modP, Model *mod_approxP);
};

struct Measures
{
    // int N;
    // double T;
    vector<vector<vector<vector<int> > > > N_AP;
    vector<vector<vector<double> > > rates;
    vector<vector<vector<vector<double> > > > rates_T;
};

struct Computation
{
    // modes
    // int p_theory, p_theory_hist;
    // vector<long int> seed_theory, seed_time;
    long int seed;
    unsigned N;
    double T;
    unsigned draw_from_theory, draw_finite_time;
    unsigned k, j;
    // int process_data;

    // parameter
//         double alpha_0, rateWnt;
    // long int seed_time;
    // int N, n_bin, border;
    // string prior;

    double draw_rate(gsl_rng *rng, Model *modP);
    double draw_sample(gsl_rng *rng, double rate, double T);
};

struct parameters_int
{
	double alpha_0, Psi_0;
	double rate_max;
	double gamma, delta;

	double gamma_approx, delta_approx;

    double I_alpha, I_beta;
};
