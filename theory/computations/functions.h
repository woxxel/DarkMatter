using namespace std;

void add_PSP(int p, double tau_I, double tau_norm, double tau_n);

void initiate_results(model *modP, simulation *simP);
void compare_approx(model *modP, model *mod_approxP);
// void find_transitions(model *modP, simulation *simP);

// int selfconsistency_f (const gsl_vector * q, void * paras, gsl_vector * f);   // dummy function for integration
double I_squared_nu(double alpha, double sigma_V, double rateWnt, double rate_max);
double I_squared_q(double alpha, double sigma_V, double q, double rate_max);
double selfcon(double alpha, double sigma_V, double rateWnt, double q, double rate_max);

double int_distribution_exact(double nu, void *params);
double int_information_distribution(double nu, void *params);
double information_fct(double nu, parameters_int *paras);

double pdf2hist(double lower, double upper, parameters paras);
// double GR_implausible_test(double lower, double upper, model mod);
// vector<double> get_density_estimate(vector<int> data, computation com, string kernel);
// void get_distribution(model *modP, simulation *simP, results *resP);
double rate_distribution(double nu, double rate_max, double gamma, double delta);

double int_shannon_entropy(double nu, void *params);
double shannon_entropy(int p, double lower, double upper, parameters paras);

double int_KL(double nu, void *params);
double KL_divergence(int p, double lower, double upper, parameters paras, parameters paras_approx);


// void bayesian_estimate_theory(model *modP, computation *comP, Results *resP);
// void bayesian_estimate_measures(measures *mesP, computation *comP, Results *resP);
// double bayes_est_prior(double nu, model *modP, string prior);
// double bayes_est_prior_measures(double nu, string prior);
//
// void draw_rates(model *modP, computation *comP, Results *resP);
// void draw_samples(computation *comP, Results *resP);
//
// vector<double> get_cdf(vector<double> p, int steps, double d_nu);
// void post_process(computation *comP, model *modP, model *mod_approxP, simulation *simP, Results *resP);

// void fill_line(int value, double *gammaP, computation com, unsigned n_iter, unsigned alpha_iter, unsigned tau_G_iter, unsigned rate_iter, unsigned p);

//! read/write stuff
// void print_state(size_t iter, gsl_multiroot_fsolver * s);

double poisson_sample(int k, double lambda, double T);
// int factorial(int n);

// unsigned max_multiple(int argc, unsigned** argv);
