#include <vector>

using namespace std;

struct parameters
{
    int Npop;
	double tau_G, tau_A, tau_N, tau_M, J;
	double eta, eps, n;
    double kappa;
    vector<double> rate, alpha_0;

    //external drive
    int drive;
    double J_0; // synaptic strength with respect to recurrent weights J
    double K_0; // average incoming number of synapses as multiple of recurrent connections number K
    double tau_0;

	vector<double> q;

	vector<double> J_E, J_I;
//         vector<double> kappa;

	vector<double> alpha_raw, alpha, sigma_V;
	vector<double> rate_max;
	vector<double> gamma, delta, chi, I_balance;
    vector<double> regions;
    vector<double> entropy, KL;
};

struct results
{
    vector<double> p_exact, p_range, p_approx;
    vector<vector<double> > rate_inf;
//         vector<double> delta, gamma, chi;
    vector<vector<vector<double> > > rate, alpha_raw, alpha, sigma_V, I_balance;

    vector<vector<vector<double> > > gamma, chi, regions;//, regions_approx;
    vector<vector<vector<double> > > q, q_approx;
    vector<vector<vector<double> > > gamma_approx, chi_approx;

    vector<unsigned> axesDim;

    double alpha_single, sigma_V_single, I_single, chi_single;
//         vector<vector<double> > q, alpha_raw, alpha, sigma_V, gamma, chi;

    vector<vector<double> > KS;
    vector<vector<vector<double> > > entropy, KL_entropy;
//         vector<double> KL_entropy;

    vector<vector<vector<double> > > rate_T;
    vector<vector<vector<int> > > N_AP;
//         vector<double> p_est_inf, p_est_T;
    vector<vector<vector<double> > > p_bayes_est;
    vector<double> p_bayes_est_measures;
    vector<double> cdf_theory;
//         vector<double> p_k;
    vector<int> p_hist;

    vector<vector<double> > trans_DM, trans_np, trans_inc, trans_imp;

    double d_nu, max_prob, factor;
    int steps;

    int axes[6] = {0,0,0,0,0,0};
};

class model
{
	public:
        parameters paras;
        void set_weights();
        void solve_selfcon(int mode_calc);
        void write_results(results * resP);
        double distribution_exact(double nu, int p);
        bool q_border1(unsigned p);
        bool q_border2(unsigned p);
        bool inconsistent();
        bool no_peak(unsigned p);
        bool implausible();

        void store_update(results *resP);

// 		model();
        void resize();

    private:

        void get_sigma_V();
        void get_alpha();
        void get_delta();
        void get_gamma();
        void get_chi();

        double nu_peak_log_full(double gamma, double delta, double rate_max);
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

        vector<double> n, alpha_0, tau_G, rateWnt, eps, eta;
        unsigned nVar;
        vector<string> order; // contains the strings of parameters in the order that the iteration should walk through

        unsigned n_iter, alpha_0_iter, tau_G_iter, rateWnt_iter, eps_iter, eta_iter;
        int mode_calc, mode_stats;
        size_t nSz, alpha_0Sz, tau_GSz, rateWntSz, epsSz, etaSz, orderSz, charSz, steps;

        vector<bool> trans_DM_found, trans_np_found;
        bool trans_imp_found, trans_inc_found;

        double y_val();

        vector<simulation_variable> vars; // initialize, such that variables are registered in order
        bool iter_status[6];
        // bool iterating = true;
        void initialize(model *modP);
        void reset_iteration();
        bool run_iteration(model *modP);
        void print_simStatus();

        void initiate_y_axis(int Npop);
        void store_results(simulation *simP, model *modP, model *mod_approxP, results *resP);
        void store_results_approx(simulation *simP, model *modP, results *resP);
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
};
