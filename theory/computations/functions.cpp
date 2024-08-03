#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include "functions.h"

// define classes

void Model_Results::initiate(unsigned sim_steps_1, unsigned sim_steps_2, unsigned sim_mode)
{
    trans_inc.resize(sim_steps_1,vector<int>(0));
    trans_imp.resize(sim_steps_1,vector<int>(0));
}

void Population_Results::initiate(unsigned sim_steps_1, unsigned sim_steps_2, unsigned sim_mode)
{
    trans_DM.resize(sim_steps_1,vector<int>(0));
    trans_np.resize(sim_steps_1,vector<int>(0));

    rate.resize(sim_steps_2,vector<double>(sim_steps_1));
    q.resize(sim_steps_2,vector<double>(sim_steps_1));
    gamma.resize(sim_steps_2,vector<double>(sim_steps_1));
    delta.resize(sim_steps_2,vector<double>(sim_steps_1));
    rate_max.resize(sim_steps_2,vector<double>(sim_steps_1));
    chi.resize(sim_steps_2,vector<double>(sim_steps_1));
    I_balance.resize(sim_steps_2,vector<double>(sim_steps_1));

    if ((sim_mode == 0) || (sim_mode == 3))
    {
        regions.resize(sim_steps_2,vector<double>(sim_steps_1));
        implausible.resize(sim_steps_2,vector<double>(sim_steps_1));
        // if (sim_mode == 3)
            // regions_approx.resize(sim_steps_2,vector<double>(sim_steps_1));
    }


    if (sim_mode == 1)
    {
        // cout << "start q size : " << resP->q.size() << endl;
        regions.resize(sim_steps_2,vector<double>(sim_steps_1));
        implausible.resize(sim_steps_2,vector<double>(sim_steps_1));
        alpha_raw.resize(sim_steps_2,vector<double>(sim_steps_1));
        alpha.resize(sim_steps_2,vector<double>(sim_steps_1));
        sigma_V.resize(sim_steps_2,vector<double>(sim_steps_1));
    }

    if ((sim_mode == 2) || (sim_mode == 3))
    {
        // q_approx.resize(sim_steps_2,vector<double>(sim_steps_1));
        // gamma_approx.resize(sim_steps_2,vector<double>(sim_steps_1));
        // chi_approx.resize(sim_steps_2,vector<double>(sim_steps_1));

        KL_entropy.resize(sim_steps_2,vector<double>(sim_steps_1));
        entropy.resize(sim_steps_2,vector<double>(sim_steps_1));

        // trans_DM_approx.resize(sim_steps_1);
        // trans_np_approx.resize(sim_steps_1);
        // trans_imp_approx.resize(sim_steps_1);
        // trans_inc_approx.resize(sim_steps_1);
    }

    if (sim_mode == 4)
    {
        // modP->infoContent.resize();
        infoContent.resize(sim_steps_2,vector<double>(sim_steps_1));
    }
}


void initiate_results(Model *modP, Simulation *simP)
{

    simP->trans_imp_found = false;
    simP->trans_inc_found = false;

    simP->trans_imp_found_approx = false;
    simP->trans_inc_found_approx = false;

    modP->results.initiate(simP->vars[0].steps,simP->vars[1].steps,simP->mode_stats);

    for (unsigned l=0; l<modP->L; l++) {
        for (unsigned p=0; p<modP->layer[l].nPop; p++)
            modP->layer[l].population[p].results.initiate(simP->vars[0].steps,simP->vars[1].steps,simP->mode_stats);
    }
}

// void get_distribution(Model *modP, Simulation *simP, results *resP)
// {
//     for (unsigned p = 0; p < modP->paras.Npop; p++)
//     {
//         resP->d_nu = modP->paras.rate_max[p]/simP->steps;
//
//         resP->p_exact[0] = 0;
//         resP->cdf_theory[0] = 0;
//         for (unsigned i=1; i<simP->steps; ++i)
//         {
//             resP->p_range[i] = i*resP->d_nu;    // should be written in some seperate function
//             resP->p_exact[i] = modP->distribution_exact(resP->p_range[i],0);
//             // resP->max_prob = max(resP->p_exact[i],resP->max_prob);
//             resP->cdf_theory[i] = resP->cdf_theory[i-1] + pdf2hist((i-1)*resP->d_nu,i*resP->d_nu,modP->paras);
//         }
//     }
// }

void compare_approx(Model *modP, Model *modP_approx)
{
    Population_Simulation *popSimP, *popSimP_approx;
    for (unsigned l = 0; l < modP->L; l++) 
    {
        for (unsigned p = 0; p < modP->layer[l].nPop; p++)
        {
            popSimP = &modP->layer[l].population[p].simulation;
            popSimP_approx = &modP_approx->layer[l].population[p].simulation;
            long double cosh_tmp = cosh(popSimP->gamma*popSimP->delta*sqrt(-2*log(pow(10,-6)/popSimP->rate_max))); // check if pdf can be calculated
            if (!isinf(cosh_tmp) && !isnan(cosh_tmp))
            {
                popSimP->KL = KL_divergence(p,0,popSimP->rate_max,popSimP,popSimP_approx);
                popSimP->entropy = shannon_entropy(p,0,popSimP->rate_max,popSimP);

            } else {
                popSimP->KL = NAN;
                popSimP->entropy = NAN;
            }
        }

    }
}

int selfconsistency_from_currents_f (const gsl_vector * vars, void * params, gsl_vector * f)
{
    Model * modP = (Model *) params;
    Population *popP;
    
    unsigned idx_pop[] = { modP->layer[0].nPop-2, modP->layer[0].nPop-1 };
    popP = &modP->layer[0].population[idx_pop[0]];

    double rate = popP->rateWnt;
    double rate_max = popP->simulation.rate_max;
    double sigma_V = popP->simulation.sigma_V;
    
    // define the two selfconsistency equations:
    long double alpha_I_sq = gsl_vector_get(vars,0);
    long double I_0 = gsl_vector_get(vars,1);
    if (isnan(I_0)) {
        return GSL_FAILURE;
    }

    
    
    long double alpha_1 = sqrt(alpha_I_sq + gsl_pow_2(modP->layer[0].population[idx_pop[0]].alpha_0));
    long double alpha_2 = sqrt(alpha_I_sq + gsl_pow_2(modP->layer[0].population[idx_pop[1]].alpha_0));

    long double selfcon_nu = 1./2. * (
        calc_first_moment(rate_max,sigma_V,alpha_1,I_0,modP->layer[0].population[idx_pop[0]].Psi_0) + 
        calc_first_moment(rate_max,sigma_V,alpha_2,I_0,modP->layer[0].population[idx_pop[1]].Psi_0)
        ) - rate;

    long double selfcon_alpha = gsl_pow_2(popP->J[idx_pop[0]])/2. * (
        calc_second_moment(rate_max,sigma_V,alpha_1,I_0,modP->layer[0].population[idx_pop[0]].Psi_0) + 
        calc_second_moment(rate_max,sigma_V,alpha_2,I_0,modP->layer[0].population[idx_pop[1]].Psi_0)
        ) - alpha_I_sq;

    // cout << "nu: (" << calc_first_moment(rate_max,sigma_V,alpha_1,I_0,modP->layer[0].population[0].Psi_0) << " , " << calc_first_moment(rate_max,sigma_V,alpha_2,I_0,modP->layer[0].population[1].Psi_0) << ")" << endl;
    
    // cout << "q: (" << calc_second_moment(rate_max,sigma_V,alpha_1,I_0,modP->layer[0].population[0].Psi_0) << " , " << calc_second_moment(rate_max,sigma_V,alpha_2,I_0,modP->layer[0].population[1].Psi_0) << ")" << endl;

    // cout << "selfcon_nu = " << selfcon_nu << endl;
    // cout << "selfcon_alpha = " << selfcon_alpha << endl;

    // gsl_vector_set(f,0,gsl_pow_2(selfcon_nu));
    // gsl_vector_set(f,1,gsl_pow_2(selfcon_alpha));
    gsl_vector_set(f,0,selfcon_nu);
    gsl_vector_set(f,1,selfcon_alpha);

    return GSL_SUCCESS;
}

long double calc_first_moment(double rate_max, double sigma_V, long double alpha, long double I_0, double Psi_0)
{
    return rate_max * sigma_V / sqrt(gsl_pow_2(alpha) + gsl_pow_2(sigma_V)) * exp(- gsl_pow_2(I_0 - Psi_0)/(2*(gsl_pow_2(alpha) + gsl_pow_2(sigma_V))));
}

long double calc_second_moment(double rate_max, double sigma_V, long double alpha, long double I_0, double Psi_0)
{
    return gsl_pow_2(rate_max) * sigma_V / sqrt(2*gsl_pow_2(alpha) + gsl_pow_2(sigma_V)) * exp(- gsl_pow_2(I_0 - Psi_0)/(2*gsl_pow_2(alpha) + gsl_pow_2(sigma_V)));
}

int selfconsistency_f_split (const gsl_vector * vars, void * params, gsl_vector * f)
{
	Model * modP = (Model *) params;

	vector<vector<double> > alpha, q_search;

    double delta_nu = gsl_vector_get(vars,0);

    unsigned p_idx = 1;
    // vector<double> q_search = {vars.begin() + 1, vars.end()};

    q_search.resize(modP->L);
    alpha.resize(modP->L);
    for (unsigned l = 0; l < modP->L; l++) {
        q_search[l].resize(modP->layer[l].nPop);
        alpha[l].resize(modP->layer[l].nPop);
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {

            modP->layer[l].population[p].simulation.rate = modP->network_rate + pow(-1,p)*delta_nu;

            q_search[l][p] = gsl_vector_get(vars,p_idx);
            p_idx++;
        }
    }

    alpha = modP->calc_alpha(q_search);
    modP->get_sigma_V();
    // cout << "alpha: 1: " << alpha[0][0] << ", 2: " << alpha[0][1] << endl;
    // cout << "sigma: 1: " << modP->layer[0].population[0].simulation.sigma_V << ", 2: " << modP->layer[0].population[1].simulation.sigma_V << endl;

    // compute and set the results
    p_idx = 0;
    for (unsigned l = 0; l < modP->L; l++) {
        gsl_vector_set (f, p_idx, selfcon_split(
            modP->network_rate + delta_nu, modP->network_rate - delta_nu, 
            modP->layer[l].population[0].Psi_0, modP->layer[l].population[1].Psi_0,
            alpha[l][0],
            modP->layer[l].population[0].simulation.sigma_V,
            modP->layer[l].population[0].simulation.rate_max)
        );
        p_idx++;
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {
            gsl_vector_set (f, p_idx, selfcon(
                alpha[l][p],
                modP->layer[l].population[p].simulation.sigma_V,
                modP->network_rate + pow(-1,p)*delta_nu,
                q_search[l][p],
                modP->layer[l].population[p].simulation.rate_max));
            p_idx++;
        }
        
    }

	return GSL_SUCCESS;
}


int selfconsistency_f (const gsl_vector * q, void * params, gsl_vector * f)
{
	Model * modP = (Model *) params;

	vector<vector<double> > alpha, q_search;

    unsigned p_idx = 0;
    q_search.resize(modP->L);
    alpha.resize(modP->L);
    for (unsigned l = 0; l < modP->L; l++) {
        q_search[l].resize(modP->layer[l].nPop);
        alpha[l].resize(modP->layer[l].nPop);
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {
            q_search[l][p] = gsl_vector_get(q,p_idx);
            p_idx++;
        }
    }

    alpha = modP->calc_alpha(q_search);

    // compute and set the results
    p_idx = 0;
    for (unsigned l = 0; l < modP->L; l++) {
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {
            gsl_vector_set (f, p_idx, selfcon(alpha[l][p],modP->layer[l].population[p].simulation.sigma_V,modP->layer[l].population[p].rateWnt,q_search[l][p],modP->layer[l].population[p].simulation.rate_max));
            p_idx++;
        }
    }

	return GSL_SUCCESS;
}


double selfcon_split(double rate_1, double rate_2, double psi_1, double psi_2, double alpha, double sigma_V, double rate_max)
{
    return sqrt(I_squared_nu(alpha,sigma_V,rate_1,rate_max)) - sqrt(I_squared_nu(alpha,sigma_V,rate_2,rate_max)) - psi_1 + psi_2;
}

double selfcon(double alpha, double sigma_V, double rate, double q, double rate_max)
{
	return I_squared_nu(alpha,sigma_V,rate,rate_max) - I_squared_q(alpha,sigma_V,q,rate_max);
}

double I_squared_nu(double alpha, double sigma_V, double rate, double rate_max)
{
    return -( gsl_pow_2(alpha) + gsl_pow_2(sigma_V) ) * log( gsl_pow_2(rate/rate_max) * (1 + gsl_pow_2(alpha / sigma_V)) );
}

double I_squared_q(double alpha, double sigma_V, double q, double rate_max)
{
    return -( gsl_pow_2(alpha) + 0.5*gsl_pow_2(sigma_V) ) * log( (gsl_pow_2(q)/gsl_pow_4(rate_max) ) * (1 + 2*gsl_pow_2(alpha / sigma_V)) );
}

// void simulation::store_results_approx(Simulation *simP, Model *modP, results * resP)
// {
// //         unsigned a = resP->gamma_approx.size() - 1;
// //
// // //         resP->rate[a].push_back(paras.rate);
//     cout << "storing" << endl;
//     for (int p = 0; p < modP->paras.Npop; p++)
//     {
//         if (mode_stats == 2)
//             resP->q_approx[p][vars[1].iter][vars[0].iter] = modP->paras.q[p];
//         cout <<
// //         cout << "gamma: " << modP->paras.gamma[0] << endl;
// //         cout << "gamma_size: " << resP->gamma_approx.size() << endl;
//         resP->gamma_approx[p][vars[1].iter][vars[0].iter] = modP->paras.gamma[p];
//         resP->chi_approx[p][vars[1].iter][vars[0].iter] = modP->paras.chi[p];
// //         resP->regions_approx[a].push_back(modP->paras.regions[0]);
//     }
//     cout << "done" << endl;
// }

// void model::store_update(results * resP)
// {
//         //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent
// }


// unsigned max_multiple(int n_args, ...)
// {
//     va_list ap
//     for(int n=1; n<argc; ++n)
//             max_tmp = argv[n] > argv[n-1] ? argv[n] : argv[n-1]
//     return max_tmp;
// }

double shannon_entropy(int p, double lower, double upper, Population_Simulation *paras)
{
//         struct parameters_int *para = (struct parameters_int *) paras;

    struct parameters_int Pparas;

    Pparas.rate_max = paras->rate_max;
    Pparas.gamma = paras->gamma;
    Pparas.delta = paras->delta;

    gsl_function p_integrate;
    p_integrate.function = &int_shannon_entropy; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
    p_integrate.params = &Pparas;
    // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
    gsl_set_error_handler_off();
    double res, err;
    gsl_integration_qag(&p_integrate, lower, upper, 1e-7, 1e-7, 10000, 6, ww, &res, &err);		//! integrate G from xth2 to infty
    gsl_integration_workspace_free (ww);

//         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
    return res;//*(upper-lower);
}

double int_shannon_entropy(double nu, void *params)
{
    struct parameters_int paras = *(struct parameters_int *) params;

    long double rho_exact = rate_distribution(nu,paras.rate_max,paras.gamma,paras.delta);

    return -rho_exact*log2(rho_exact);
}


double KL_divergence(int p, double lower, double upper, Population_Simulation *paras, Population_Simulation *paras_approx)
{
//         struct parameters_int *para = (struct parameters_int *) paras;

    struct parameters_int Pparas;

    Pparas.rate_max = paras->rate_max;
    Pparas.gamma = paras->gamma;
    Pparas.delta = paras->delta;

    Pparas.gamma_approx = paras_approx->gamma;
    Pparas.delta_approx = paras_approx->delta;

    gsl_function p_integrate;
    p_integrate.function = &int_KL; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
    p_integrate.params = &Pparas;
    // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
    gsl_set_error_handler_off();
    double res, err;
    gsl_integration_qag(&p_integrate, lower, upper, 1e-7, 1e-7, 10000, 6, ww, &res, &err);		//! integrate G from xth2 to infty
    gsl_integration_workspace_free (ww);

//         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
    return res;//*(upper-lower);
}

double int_KL(double nu, void *params)
{
        struct parameters_int paras = *(struct parameters_int *) params;
    //         struct parameters_int *paras = (struct parameters_int *) parameters;

        long double rho_exact = rate_distribution(nu,paras.rate_max,paras.gamma,paras.delta);
        long double rho_approx = rate_distribution(nu,paras.rate_max,paras.gamma_approx,paras.delta_approx);

        return rho_exact*log2(rho_exact/rho_approx);
}

double rate_distribution(double nu, double rate_max, double gamma, double delta)
{
    double rate_ratio = nu/rate_max;

    // cout << "in rate distribution" << endl;
    // cout << "distr paras: nu=" << nu << ", rate max=" << rate_max << ", rate_ratio=" << rate_ratio << ", gamma="<< gamma << ", delta=" << delta << endl;

    if (nu == 0)
    {
        return gamma > 1 ? 0 : 1.0/0.0;
    } else {
        return gamma/(rate_max*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(delta)/2)*pow(rate_ratio,gsl_pow_2(gamma)-1)*cosh(gamma*delta*sqrt(-2*log(rate_ratio)));
    }
}


double pdf2hist(double lower, double upper, parameters paras)
{
//         struct parameters_int *para = (struct parameters_int *) paras;

    struct parameters_int Pparas;

    Pparas.rate_max = paras.rate_max[0];
    Pparas.gamma = paras.gamma[0];
    Pparas.delta = paras.delta[0];

    gsl_function p_integrate;
    p_integrate.function = &int_distribution_exact; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
    p_integrate.params = &Pparas;
    // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
    double res, err;
    gsl_integration_qag(&p_integrate, lower, upper, 1e-7, 1e-7, 10000, 6, ww, &res, &err);
    gsl_integration_workspace_free (ww);

//         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
    return res*(upper-lower);
}

double int_distribution_exact(double nu, void *params)
{
    struct parameters_int paras = *(struct parameters_int *) params;
    return rate_distribution(nu,paras.rate_max,paras.gamma,paras.delta);
}

double information_fct(double nu, parameters_int *paras)
{
    // cout << "paras in info: nu=" << nu << ", nu0=" << nu0 << ", zeta=" << zeta << ", c=" << c << endl;
    // cout << "evaluating info: " << c*pow(nu - nu0,zeta) << endl;

    // if (nu > 0.9*nu0) return 0;

    // double base = 0.5;
    // double a = zeta>0 ? zeta : base;
    // double b = zeta<0 ? -zeta : base;
    double x = nu/paras->rate_max;

    // cout << "paras in info: x=" << x << ", a=" << a << ", b=" << b << ", beta(a,b) = " << std::beta(a,b) << ", beta_distr = " << (pow(x,a-1) * pow(1-x,b-1) / std::beta(a,b)) << endl;


    return pow(x,paras->I_alpha-1) * pow(1-x,paras->I_beta-1) / std::beta(paras->I_alpha,paras->I_beta);

    // return c*pow(nu - nu0,zeta);
}

double int_information_distribution(double nu, void *params)
{
    struct parameters_int paras = *(struct parameters_int *) params;

    // cout << "nu: " << nu/paras.rate_max << ", p(nu)=" << rate_distribution(nu,paras.rate_max,paras.gamma,paras.delta) << "; I(nu)=" << information_fct(nu,paras.rate_max,paras.zeta) << endl;

    return information_fct(nu,&paras) * nu * rate_distribution(nu,paras.rate_max,paras.gamma,paras.delta);
}
