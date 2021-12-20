#include <vector>
#include <iostream>

using namespace std;

class simulation;
class model;

struct results
{
    vector<vector<vector<double> > > rate, q, q_approx, alpha_raw, alpha, sigma_V, gamma, gamma_approx, chi, chi_approx, delta, I_balance, regions, regions_approx, entropy, KL_entropy;

    vector<vector<vector<double> > > infoContent;

    vector<unsigned> axesDim;

    vector<vector<double> > KS;
    // vector<double> KL_entropy;

    vector<vector<vector<double> > > rate_T;
    vector<vector<vector<int> > > N_AP;
    vector<vector<double> > rate_inf;
    // vector<double> p_est_inf, p_est_T;
    vector<vector<vector<double> > > p_bayes_est;
    vector<double> p_bayes_est_measures;
    vector<double> cdf_theory;
    // vector<double> p_k;
    vector<int> p_hist;

    vector<vector<double> > trans_DM, trans_np, trans_inc, trans_imp;
    vector<vector<double> > trans_DM_approx, trans_np_approx, trans_inc_approx, trans_imp_approx;

    double d_nu, max_prob, factor;
    int steps;

    int axes[7] = {0,0,0,0,0,0,0};
};


struct parameters
{
    unsigned Npop;
	double tau_G, tau_A, tau_N, tau_M, J;
	double eta, eps, n;
    double zeta;
    double kappa;
    vector<double> rate, q, alpha_0;

    //external drive
    int drive;
    double J_0; // synaptic strength with respect to recurrent weights J
    double K_0; // average incoming number of synapses as multiple of recurrent connections number K
    double tau_0;

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

    unsigned nZeta;
    double minZeta, maxZeta;
};

class model
{
	public:
        parameters paras;

        vector<double> trans_DM, trans_np, trans_inc, trans_imp;
        vector<bool> trans_DM_found, trans_np_found;
        bool trans_imp_found, trans_inc_found;

        vector<double> infoContent;

        vector<bool> in_DM, in_np;
        bool in_imp, in_inc;

        void set_weights();
        void solve_selfcon(int mode_calc);
        void write_results(results * resP);
        double distribution_exact(double nu, int p);
        bool q_border1(unsigned p);
        bool q_border2(unsigned p);

        bool DM_border(unsigned p);
        bool inconsistent();
        bool no_peak(unsigned p);
        bool implausible();

        void integrate_information(info_paras infoParas);

        void find_transitions(simulation *simP, results *resP);
        void store_update(results *resP);


// 		model();
        void resize();

    private:

        void get_sigma_V();
        void get_alpha();
        void get_delta();
        void get_gamma();
        void get_chi();

        double nu_peak_log_full(unsigned p);
};

struct simulation_variable
{
    double *valP, *paraP;

    unsigned paraSz;

    string name;
    unsigned iter, steps;

    // void (model::*onUpdate)();    // function which is called every time the variable updates

    void initialize(double *modVarP, unsigned modVarSz, double *varP, unsigned varSz, string varName);
    bool iterate();
    void print_status();

};


class simulation
{
    public:

        info_paras infoParas;

        vector<double> n, alpha_0, tau_G, rateWnt, eps, eta, zeta;
        unsigned nVar;
        vector<string> order; // contains the strings of parameters in the order that the iteration should walk through

        // unsigned n_iter, alpha_0_iter, tau_G_iter, rateWnt_iter, eps_iter, eta_iter;
        int mode_calc, mode_stats;
        size_t nSz, alpha_0Sz, tau_GSz, rateWntSz, epsSz, etaSz, zetaSz, orderSz, charSz, steps;

        vector<bool> trans_DM_found, trans_np_found, trans_DM_found_approx, trans_np_found_approx;
        bool trans_imp_found, trans_inc_found, trans_imp_found_approx, trans_inc_found_approx;

        double y_val();

        vector<simulation_variable> vars; // initialize, such that variables are registered in order
        bool iter_status[6];
        // bool iterating = true;
        void initialize(model *modP);
        void reset_iteration();
        bool run_iteration(model *modP);
        void print_simStatus();

        void initiate_y_axis(model *modP);
        void store_results(model *modP, model *mod_approxP, results *resP);
};

struct computation
{
    // modes
    int p_theory, p_theory_hist;
    int draw_from_theory, draw_finite_time;
    int process_data;

    // parameter
//         double alpha_0, rateWnt;
    vector<long int> seed_theory, seed_time;
    int N, n_bin, border;
    double T;
    string prior;
    int k, j;
};

struct measures
{
    int N;
    vector<int> N_AP;
    vector<double> rates;
    double T;
};

struct parameters_int
{
	double alpha_0;
	double rate_max;
	double gamma, delta;

	double gamma_approx, delta_approx;

    double zeta;
};
