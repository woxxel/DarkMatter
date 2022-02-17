// #include <math.h>
// #include <stdarg.h>
// #include <gsl/gsl_randist.h>
// #include "gsl/gsl_rng.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include "functions.h"

// define classes
void initiate_results(model *modP, simulation *simP, results *resP)
{
    simP->trans_DM_found.resize(modP->paras.Npop,false);
    simP->trans_np_found.resize(modP->paras.Npop,false);
    simP->trans_imp_found = false;
    simP->trans_inc_found = false;

    simP->trans_DM_found_approx.resize(modP->paras.Npop,false);
    simP->trans_np_found_approx.resize(modP->paras.Npop,false);
    simP->trans_imp_found_approx = false;
    simP->trans_inc_found_approx = false;

    // now: size of border vectors
    resP->trans_DM.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_np.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_inc.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_imp.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));


    resP->rate.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->q.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->gamma.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->chi.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->delta.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->I_balance.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

    if ((simP->mode_stats == 0) || (simP->mode_stats == 3))
    {
        resP->regions.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        if (simP->mode_stats == 3)
            resP->regions_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    }


    if (simP->mode_stats == 1)
    {
        // cout << "start q size : " << resP->q.size() << endl;
        resP->alpha_raw.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->alpha.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->sigma_V.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    }

    if ((simP->mode_stats == 2) || (simP->mode_stats == 3))
    {
        resP->q_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->steps,vector<double>(simP->steps)));
        resP->gamma_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->chi_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

        resP->KL_entropy.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->entropy.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

        resP->trans_DM_approx.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
        resP->trans_np_approx.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
        resP->trans_imp_approx.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
        resP->trans_inc_approx.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    }

    if (simP->mode_stats == 4)
    {
        modP->infoContent.resize(modP->paras.Npop);
        resP->infoContent.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    }
}

// void get_distribution(model *modP, simulation *simP, results *resP)
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

void compare_approx(model *modP, model *mod_approxP)
{
    for (unsigned p = 0; p < modP->paras.Npop; p++)
    {
        long double cosh_tmp = cosh(modP->paras.gamma[p]*modP->paras.delta[p]*sqrt(-2*log(pow(10,-6)/modP->paras.rate_max[p]))); // check if pdf can be calculated
        if (!isinf(cosh_tmp) && !isnan(cosh_tmp))
        {
            modP->paras.KL[p] = KL_divergence(p,0,modP->paras.rate_max[p],modP->paras,mod_approxP->paras);
            modP->paras.entropy[p] = shannon_entropy(p,0,modP->paras.rate_max[p],modP->paras);

        } else {
            modP->paras.KL[p] = NAN;
            modP->paras.entropy[p] = NAN;
        }
    }
}


int selfconsistency_f (const gsl_vector * q, void * paras, gsl_vector * f)
{
	parameters * paraP = (parameters *) paras;

	vector<double> q_search(paraP->Npop);

	for (unsigned p = 0; p < paraP->Npop; p++)
        q_search[p] = gsl_vector_get(q,p);

	for (unsigned p = 0; p < paraP->Npop; p++)
	{
		// get alpha
        double alpha_sq = gsl_pow_2(paraP->J_I[p]) * q_search[0];

        if (paraP->Npop > 1)
            alpha_sq += gsl_pow_2(paraP->J_E[p]) * paraP->kappa * q_search[1];

        double alpha_sq_0;
        if (paraP->drive == 2)
        {
            // quenched variance from afferent, spiking drive (gauss distributed synapse numbers)
            alpha_sq_0 = sqrt(1./paraP->K_0) * gsl_pow_2(paraP->J_I[p]) * gsl_pow_2(paraP->rate[p]);
        }
        else
            alpha_sq_0 = 0;

        paraP->alpha[p] = sqrt( alpha_sq + alpha_sq_0 + gsl_pow_2(paraP->alpha_0[p]) );
		// set the function
		gsl_vector_set (f, p, selfcon(paraP->alpha[p],paraP->sigma_V[p],paraP->rate[p],q_search[p],paraP->rate_max[p]));
	}
	return GSL_SUCCESS;
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

// void simulation::store_results_approx(simulation *simP, model *modP, results * resP)
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

double shannon_entropy(int p, double lower, double upper, parameters paras)
{
//         struct parameters_int *para = (struct parameters_int *) paras;

    struct parameters_int Pparas;

    Pparas.rate_max = paras.rate_max[p];
    Pparas.gamma = paras.gamma[p];
    Pparas.delta = paras.delta[p];

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


double KL_divergence(int p, double lower, double upper, parameters paras, parameters paras_approx)
{
//         struct parameters_int *para = (struct parameters_int *) paras;

    struct parameters_int Pparas;

    Pparas.rate_max = paras.rate_max[p];
    Pparas.gamma = paras.gamma[p];
    Pparas.delta = paras.delta[p];

    Pparas.gamma_approx = paras_approx.gamma[p];
    Pparas.delta_approx = paras_approx.delta[p];

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
