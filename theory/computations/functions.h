using namespace std;

void add_PSP(int p, double tau_I, double tau_norm, double tau_n);

void initiate_results(Model *modP, Simulation *simP);
void compare_approx(Model *modP, Model *mod_approxP);
// void find_transitions(Model *modP, Simulation *simP);

// int selfconsistency_f (const gsl_vector * q, void * paras, gsl_vector * f);   // dummy function for integration
double I_squared_nu(double alpha, double sigma_V, double rateWnt, double rate_max);
double I_squared_q(double alpha, double sigma_V, double q, double rate_max);
double selfcon(double alpha, double sigma_V, double rateWnt, double q, double rate_max);
double selfcon_split(double rate_1, double rate_2, double psi_1, double psi_2, double alpha, double sigma_V, double rate_max);

long double calc_first_moment(double rate_max, double sigma_V, long double alpha, long double I_0, double Psi_0);
long double calc_second_moment(double rate_max, double sigma_V, long double alpha, long double I_0, double Psi_0);


double int_distribution_exact(double nu, void *params);
double int_information_distribution(double nu, void *params);
double information_fct(double nu, parameters_int *paras);

double pdf2hist(double lower, double upper, parameters paras);
// double GR_implausible_test(double lower, double upper, model mod);
// vector<double> get_density_estimate(vector<int> data, computation com, string kernel);
// void get_distribution(Model *modP, Simulation *simP, results *resP);
double rate_distribution(double nu, double rate_max, double gamma, double delta);

double int_shannon_entropy(double nu, void *params);
double shannon_entropy(int p, double lower, double upper, Population_Simulation *paras);

double int_KL(double nu, void *params);
double KL_divergence(int p, double lower, double upper, Population_Simulation *paras, Population_Simulation *paras_approx);


// void bayesian_estimate_theory(Model *modP, computation *comP, Results *resP);
// void bayesian_estimate_measures(measures *mesP, computation *comP, Results *resP);
// double bayes_est_prior(double nu, Model *modP, string prior);
// double bayes_est_prior_measures(double nu, string prior);
//
// void draw_rates(Model *modP, computation *comP, Results *resP);
// void draw_samples(computation *comP, Results *resP);
//
// vector<double> get_cdf(vector<double> p, int steps, double d_nu);
// void post_process(computation *comP, Model *modP, model *mod_approxP, Simulation *simP, Results *resP);

// void fill_line(int value, double *gammaP, computation com, unsigned n_iter, unsigned alpha_iter, unsigned tau_G_iter, unsigned rate_iter, unsigned p);

//! read/write stuff
// void print_state(size_t iter, gsl_multiroot_fsolver * s);

double poisson_sample(int k, double lambda, double T);
// int factorial(int n);

// unsigned max_multiple(int argc, unsigned** argv);
