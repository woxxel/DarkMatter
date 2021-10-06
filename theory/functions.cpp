#include <vector>
#include <math.h>
#include <stdarg.h>
// #include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_randist.h>
#include "gsl/gsl_rng.h"
#include <gsl/gsl_integration.h>

#include "functions.h"

// define classes

void initiate_results(model *modP, simulation *simP, results *resP)
{
    // first of: steps!
    unsigned axes_ct = 2;
    resP->axesDim.resize(6);

    resP->axesDim[0] = simP->alpha_0Sz;
    if (simP->alpha_0Sz >= 10)
    {
        resP->axes[0] = axes_ct;
        axes_ct--;
    }
    resP->axesDim[1] = simP->nSz;
    if (simP->nSz >= 10)
    {
        resP->axes[1] = axes_ct;
        axes_ct--;
    }
    resP->axesDim[2] = simP->epsSz;
    if (simP->epsSz >= 10)
    {
        resP->axes[2] = axes_ct;
        axes_ct--;
    }
    resP->axesDim[3] = simP->etaSz;
    if (simP->etaSz >= 10)
    {
        resP->axes[3] = axes_ct;
        axes_ct--;
    }
    resP->axesDim[4] = simP->tau_GSz;
    if (simP->tau_GSz >= 10)
    {
        resP->axes[4] = axes_ct;
        if (axes_ct > 0) axes_ct--;
    }
    resP->axesDim[5] = simP->rateWntSz;
    if (simP->rateWntSz >= 10)
    {
        resP->axes[5] = axes_ct;
        if (axes_ct > 0) axes_ct--;
    }

    for (unsigned i=0;i<resP->axesDim.size();++i)
    {
//                 cout << "axes '" << i << "': " << resP->axes[i] << endl;
        for (int a=0; a<2; ++a)
            if (resP->axesDim[i] > simP->max_ax[a])
            {
                if (a == 0)
                    simP->max_ax[1] = simP->max_ax[0];

                simP->max_ax[a] = resP->axesDim[i];
//                         cout << "new max " << a << ": " << simP->max_ax[a] << endl;
                break;
            }
    }
    simP->steps = simP->max_ax[0];
//         for (int a=0; a<2; ++a)
//             cout << "max " << a << ": " << simP->max_ax[a] << endl;

//         nth_element(resP->axesDim.begin(),resP->axesDim.begin()+1,resP->axesDim.end());
//         cout << "scnd largst element: " << resP->axesDim[1] << endl;

    if (axes_ct == 2)
    {
        simP->initiate_y_axis(modP);
        resP->steps = 1000;
    }
    else
        resP->steps = simP->steps;

//         cout << "mode_stats = " << simP->mode_stats << endl;
//         cout << "Npop = " << modP->paras.Npop << endl;


    if (axes_ct == 0)
    {
        simP->trans_DM_found.resize(modP->paras.Npop,false);
        simP->trans_np_found.resize(modP->paras.Npop,false);
        simP->trans_imp_found = false;
        simP->trans_inc_found = false;

        // now: size of border vectors
        resP->trans_DM.resize(modP->paras.Npop,vector<double>(resP->steps));
        resP->trans_np.resize(modP->paras.Npop,vector<double>(resP->steps));
        resP->trans_inc.resize(modP->paras.Npop,vector<double>(resP->steps));
        resP->trans_imp.resize(modP->paras.Npop,vector<double>(resP->steps));

        resP->rate.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));

        if (simP->mode_stats != 2)
        {
            resP->gamma.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->chi.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
        }


        if ((simP->mode_stats == 0) || (simP->mode_stats == 4))
            resP->regions.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));



        if (simP->mode_stats == 4)
        {
            resP->gamma_approx.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->chi_approx.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
//                 resP->regions_approx.push_back(vector<double>());

            resP->KL_entropy.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->entropy.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
        }
    }
    else if (axes_ct == 1)
    {
        simP->x_iter = -1;
        simP->y_iter = 0;
        if (simP->mode_stats == 1)
        {
            int dim1 = simP->max_ax[1];

            simP->trans_imp_found = false;
            resP->trans_imp.resize(modP->paras.Npop,vector<double>(resP->steps));

            resP->rate.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
            resP->gamma.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
            resP->chi.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));

            resP->q.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
//                         cout << "start q size : " << resP->q.size() << endl;
            resP->alpha_raw.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
            resP->alpha.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
            resP->sigma_V.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
            resP->I_balance.resize(modP->paras.Npop,vector<vector<double> >(dim1,vector<double>(resP->steps)));
        }
    }
    else
    {
        if (simP->mode_stats == 3)
        {
            resP->rate.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->q.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->q_approx.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));

            resP->gamma.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->gamma_approx.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->chi.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->chi_approx.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));

            resP->KL_entropy.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
            resP->entropy.resize(modP->paras.Npop,vector<vector<double> >(resP->steps,vector<double>(resP->steps)));
        }
    }


    if ((simP->mode_stats == 3) || (simP->mode_stats == 4))
    {
        resP->p_range.resize(resP->steps);
        resP->p_exact.resize(resP->steps);
        resP->cdf_theory.resize(resP->steps);
    }
}


void model::resize()
{
    // resize all to having p populations
    paras.rate.resize(paras.Npop);
    paras.alpha_0.resize(paras.Npop);

	paras.q.resize(paras.Npop);

	paras.J_E.resize(paras.Npop);
	paras.J_I.resize(paras.Npop);
//         paras.kappa.resize(paras.Npop);

    paras.alpha_raw.resize(paras.Npop);
	paras.alpha.resize(paras.Npop);
	paras.sigma_V.resize(paras.Npop);
	paras.rate_max.resize(paras.Npop);

	paras.gamma.resize(paras.Npop);
	paras.delta.resize(paras.Npop);
    paras.I_balance.resize(paras.Npop);
	paras.chi.resize(paras.Npop);
    paras.regions.resize(paras.Npop);

    paras.KL.resize(paras.Npop);
    paras.entropy.resize(paras.Npop);

}

void model::set_weights()
{
//         cout << "number of populations: " << paras.Npop << endl;

        if (paras.drive == 2)
                paras.J_0 = paras.tau_M;

        if (paras.Npop == 1)
        {
                paras.J_I[0] = paras.tau_M;				// J_II
                paras.J_I[1] = paras.tau_M;		// J_EI

                paras.J_E[0] = 0;			                                // J_IE
                paras.J_E[1] = 0;					// J_EE

        }
        else if (paras.Npop == 2)
        {
                // watch out for indexing:
                // inhibitory population: index 0
                // excitatory population: index 1
                paras.J_I[0] = sqrt(1 - gsl_pow_2(paras.eps)) * paras.tau_M;				// J_II
                paras.J_I[1] = sqrt(1 - gsl_pow_2(paras.eta * paras.eps)) * paras.tau_M;		// J_EI

                paras.J_E[0] = paras.eps * paras.tau_M;			                                // J_IE
                paras.J_E[1] = paras.eta * paras.eps * paras.tau_M;					// J_EE

                //! for p=1, the excitatory population receives inhibition, but is decoupled, such that it gives no feedback
//                 cout << "J_I: " << paras.J_I[0] << "," << paras.J_I[1] << endl;
//                 cout << "J_E: " << paras.J_E[0] << "," << paras.J_E[1] << endl;
        }
        else
        {
                cout << "There is no algorithm set yet to assign weights for more than 2 populations. Please fix that before you continue!" << endl;
                exit(0);
        }
}

void model::get_sigma_V()
{
	// iterate over both populations
	for (int p = 0; p < paras.Npop; p++)
	{
		// get sigma_V
//                 cout << "tau_G: " << paras.tau_G << ", tau_A: " << paras.tau_A << ", tau_N: " << paras.tau_N << ", tau_M: " << paras.tau_M << endl;
                // from excitatory AMPA synapses
                double var_V_A = gsl_pow_2(paras.J_E[p]) * paras.kappa * paras.rate[p] / (paras.tau_A + paras.tau_M) * ( gsl_pow_2(1-paras.n)/2 + (1-paras.n)*paras.n*paras.tau_A / (paras.tau_A + paras.tau_N) );
                // from excitatory NMDA synapses
		double var_V_N = gsl_pow_2(paras.J_E[p]) * paras.kappa * paras.rate[p] / (paras.tau_N + paras.tau_M) * ( paras.n*paras.n/2 + (1-paras.n)*paras.n*paras.tau_N / (paras.tau_A + paras.tau_N) );
                // from inhibitory GABA synapses
		double var_V_G = gsl_pow_2(paras.J_I[p]) * paras.rate[p] * 0.5 / (paras.tau_G + paras.tau_M);

                // from external drive
//                 cout << "J_I = " << paras.J_I[p] << ", J_0 = " << paras.J_0 << endl;

                // total
		double var_V = var_V_A + var_V_N + var_V_G;
//                 cout << "var total: " << var_V << ", from external sources: " << var_V_0 << endl;

		double var_V_dot = var_V_A / (paras.tau_A * paras.tau_M) + var_V_N / (paras.tau_N * paras.tau_M) + var_V_G / (paras.tau_G * paras.tau_M);

                if (paras.drive == 2)
                {
                        double var_V_0 = sqrt(1/paras.K_0) * paras.J_0 * paras.J_I[p] * paras.rate[p] * 0.5 / (paras.tau_0 + paras.tau_M);
                        var_V += var_V_0;
                        var_V_dot += var_V_0 / (paras.tau_0 * paras.tau_M);
                }
		paras.sigma_V[p] = sqrt(var_V);

		// and the maximum firing rate response
		paras.rate_max[p] = sqrt(var_V_dot / var_V) / (2 * M_PI);
//                 cout << "sigma population " << p << ": " << paras.sigma_V[p] << endl;

	}
}

void model::get_alpha()
{
        for (int p = 0; p < paras.Npop; p++)
        {
                double alpha_sq = gsl_pow_2(paras.J_I[p]) * paras.q[0];

                if (paras.Npop > 1)
                        alpha_sq += gsl_pow_2(paras.J_E[p]) * paras.kappa * paras.q[1];

                double alpha_sq_0;
                if (paras.drive == 2)
                {
                        // quenched variance from afferent, spiking drive (gauss distributed synapse numbers)
                        // after substitution of v_0 = nu_I * J/J_0 * sqrt(K/K_0)
                        alpha_sq_0 = 1/paras.K_0 * gsl_pow_2(paras.J_I[p]) * gsl_pow_2(paras.rate[p]);
                }
                else
                        alpha_sq_0 = 0;

                paras.alpha_raw[p] = sqrt( alpha_sq + alpha_sq_0);
                paras.alpha[p] = sqrt( alpha_sq + alpha_sq_0 + gsl_pow_2(paras.alpha_0[p]));
//                 cout << "alpha total: " << gsl_pow_2(paras.alpha[p]) << ", from external sources: " << alpha_sq_0 << endl;
        }
}

void model::get_delta()
{
        for (int p = 0; p < paras.Npop; p++)
        {
                double I_balance = sqrt( I_squared_nu(paras.alpha[p],paras.sigma_V[p],paras.rate[p],paras.rate_max[p]) );
                paras.I_balance[p] = I_balance;
                paras.delta[p] = I_balance/paras.alpha[p];
        }
}

//! all those could have an added boolean, showing whether they were calculated yet, and only evaluate, if not. however, this would require boolean to be updated, whenever some variable is changed
void model::get_gamma()
{
        for (int p = 0; p < paras.Npop; p++)
                paras.gamma[p] = paras.sigma_V[p]/paras.alpha[p];
}

void model::get_chi()
{
        for (int p = 0; p < paras.Npop; p++)
	{
                double nu_peak_log_I = nu_peak_log_full(paras.gamma[p],paras.delta[p],paras.rate_max[p]);
                paras.chi[p] = -log10(exp(1)) * nu_peak_log_I + log10(paras.rate[p]);
        }
}

double model::nu_peak_log_full(double gamma, double delta, double rate_max)
{
        return log(rate_max) - (gsl_pow_2(gamma*delta)-2*(gsl_pow_2(gamma)- 1) + gamma*delta*sqrt(gsl_pow_2(gamma*delta)-4*(gsl_pow_2(gamma)-1))) /(4*gsl_pow_2(gsl_pow_2(gamma) - 1));
}

double model::distribution_exact(double nu, int p)
{
        //! should not be evaluated, if gamma*delta > 200 or so, as cosh -> infty
        double rate_ratio = nu/paras.rate_max[p];
//         cout << "ratio: " << rate_ratio << endl;
//         cout << "log: " << log(rate_ratio) << endl;
//         cout << "cosh: " << cosh(paras.gamma[p]*paras.delta[p]*sqrt(-2*log(rate_ratio))) << endl;
        return paras.gamma[p]/(paras.rate_max[p]*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(paras.delta[p])/2)*pow(rate_ratio,gsl_pow_2(paras.gamma[p])-1)*cosh(paras.gamma[p]*paras.delta[p]*sqrt(-2*log(rate_ratio)));
}

void model::solve_selfcon(int mode_calc)
{
        // initiate variables
//         set_weights();                             // set parameters

        get_sigma_V();

//         cout << "solving the selfconsistency equations with the " << mode_calc << " mode." << endl;

        //! solving selfconsistency equations(s)
        if (mode_calc == 0) // exact solution
        {
                // numeric root finding to get exact solution
                // set solver parameter
                gsl_vector *q_guess = gsl_vector_alloc (paras.Npop);
                gsl_multiroot_function F;

//                 int f_tmp = this->selfconsistency_f;
//                 F.f = &f_tmp;
                F.f = &selfconsistency_f;
                F.n = paras.Npop;
                F.params = &paras;
                for (int p = 0; p < paras.Npop; p++)
                {
                        gsl_vector_set(q_guess, p, gsl_pow_2(paras.rate[p]));	         // first guess
//                         cout << "q guess (p=" << p << "): " << gsl_pow_2(paras.rate[p]) << endl;
                }
                // initiate the solver (dnewton)
                int status;
                const gsl_multiroot_fsolver_type *Tsolv = gsl_multiroot_fsolver_dnewton;
                gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (Tsolv, paras.Npop);
                gsl_multiroot_fsolver_set (s, &F, q_guess);
                // iterate the solver
                size_t iter = 0;
                do
                {
                        iter++;
                        status = gsl_multiroot_fsolver_iterate (s);

                        if (status) break;	// check if solver is stuck

                        status = gsl_multiroot_test_residual (s->f, 1e-7);
                }
                while (status == GSL_CONTINUE && iter < 100);

                //! post processing
                for (int p = 0; p < paras.Npop; p++)
                {
                    paras.q[p] = gsl_vector_get(s->x,p);
//                     cout << "q(" << p << ",exact)=" << paras.q[p] << endl;
                }
                gsl_multiroot_fsolver_free(s);
        }
        if (mode_calc == 1) // analytic approximation
        {
//                 cout << "apply approximation" << endl;
                // apply "xyz"-approximation to obtain selfconsistent second moment
                for (int p = 0; p < paras.Npop; ++p)
                {
                    double tau_q = 2*(paras.tau_M + paras.tau_G);
//                     cout << "tau_M=" << paras.tau_M << ", tau_G=" << paras.tau_G << endl;
                    double potent = tau_q*(gsl_pow_2(paras.rate[p]) + paras.alpha_0[p]/gsl_pow_2(paras.J_I[p]))/(paras.rate[p] + 2*tau_q*(gsl_pow_2(paras.rate[p]) + paras.alpha_0[p]/gsl_pow_2(paras.J_I[p])));


                    paras.q[p] = 1./sqrt(1+2*tau_q*(paras.rate[p] + paras.alpha_0[p]/(gsl_pow_2(paras.J_I[p])*paras.rate[p])))*pow(1+tau_q*(paras.rate[p] + paras.alpha_0[p]/(gsl_pow_2(paras.J_I[p])*paras.rate[p])),1-potent)*pow(paras.rate_max[p],2*potent)*pow(paras.rate[p],2-2*potent);
//                     cout << "q(" << p << ",approx)=" << paras.q[p] << endl;
                }
        }

        get_alpha();
        get_delta();
        get_gamma();
        get_chi();

//         cout << "rate = " << paras.rate[0] << ":\t q=" << paras.q[0] << "\t alpha=" << paras.alpha[0] << " ,\t sigma=" << paras.sigma_V[0] << "delta=" << paras.delta[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
}

bool model::q_border1(unsigned p)
{
        return (1 < paras.rate[p]/paras.rate_max[p] * sqrt(gsl_pow_2(paras.alpha[p]) + gsl_pow_2(paras.sigma_V[p]))/paras.sigma_V[p]);
}

bool model::q_border2(unsigned p)
{
        return (1 < paras.q[p]/gsl_pow_2(paras.rate_max[p]) * sqrt(2 * gsl_pow_2(paras.alpha[p]) + gsl_pow_2(paras.sigma_V[p]))/paras.sigma_V[p]);
}

bool model::inconsistent()
{
        for (int p = 0; p < paras.Npop; p++)
                if (q_border1(p) || q_border2(p))
                        return true;
        return false;
}

bool model::no_peak(unsigned p)
{

        return (gsl_pow_2(paras.gamma[p] * paras.delta[p])-4*(gsl_pow_2(paras.gamma[p]) - 1) < 0);
}

bool model::implausible()
{
        double lower = 0.9;
        double upper = 1;

//         cout << "rate max: " << paras.rate_max[0] << ", gamma: " << paras.gamma[0] << ", delta: " << paras.delta[0] << endl;
        for (int p = 0; p < paras.Npop; p++)
        {
                struct parameters_int Pparas;

                Pparas.rate_max = paras.rate_max[p];
                Pparas.gamma = paras.gamma[p];
                Pparas.delta = paras.delta[p];

                gsl_function p_integrate;
                p_integrate.function = &int_distribution_exact; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
                p_integrate.params = &Pparas;
                // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
                gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
                gsl_set_error_handler_off();
                double res, err;
        //         try
        //         {
                gsl_integration_qags(&p_integrate, lower*paras.rate_max[p], upper*paras.rate_max[p], 1e-7, 1e-7, 10000, ww, &res, &err);		//! integrate G from xth2 to infty
        //         } catch (const std::exception& e){
        //             cout << "whatup?" << endl;
        //         } catch (...)
        //         {
        //             cout << "error while integrating..." << endl;
        //         }
                gsl_integration_workspace_free (ww);

        //         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
                if (res > 0.1)
                        return true;
        }
        return false;
}



int selfconsistency_f (const gsl_vector * q, void * paras, gsl_vector * f)
{
	parameters * paraP = (parameters *) paras;

	vector<double> q_search(paraP->Npop);

	for (int p = 0; p < paraP->Npop; p++)
                q_search[p] = gsl_vector_get(q,p);

	for (int p = 0; p < paraP->Npop; p++)
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



double simulation::y_val()
{
        if (tau_GSz > 1)
            return tau_G[tau_G_iter];
        else if (alpha_0Sz > 1)
            return alpha_0[alpha_0_iter];
        else if (nSz > 1)
            return n[n_iter];
        else
            return 0;
//         {
//             cout << "no y-axis could be found - quitting..." << endl;
//             exit(0);
//         }
}

void simulation::initiate_y_axis(model *modP)
{
//         cout << "initiating y-axis..." << endl;

    if (mode_stats == 0 || mode_stats == 2 || mode_stats == 4)
    {
        x_iter = -1;
        y_iter++; // keep track of current y-step

        for (int p = 0; p < modP->paras.Npop; p++)
        {
            trans_DM_found[p] = false;
            trans_np_found[p] = false;
        }
    }
    else if (mode_stats == 3)// || mode_stats == 1)
        y_iter = 0;

    trans_inc_found = false;
    trans_imp_found = false;
}

void simulation::store_results(simulation *simP, model * modP, results * resP)
{
//         unsigned a = resP->rate.size() - 1;
//         cout << "size = " << a << endl;
//         unsigned size = resP->rate[a].size()+1;
//         cout << "rate=" << paras.rate << " ,\t q=" << paras.q[0] << " ,\t alpha=" << paras.alpha[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
        //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent

        for (int p = 0; p < modP->paras.Npop; p++)
        {
//                 cout << "in" << endl;
//                 cout << "rate: " << modP->paras.rate[p] << endl;
//                 cout << " x/y : " << x_iter << "/" << y_iter << endl;
//                 cout << "size rate: " << resP->rate.size() << endl;
                resP->rate[p][y_iter][x_iter] = modP->paras.rate[p];

//                 cout << " x/y : " << x_iter << "/" << y_iter << endl;
//                 cout << "rate: " << modP->paras.rate[p] << endl;
//                 cout << "storing gamma / chi" << endl;
                if (mode_stats != 2)
                {
//                         cout << "gamma: " << modP->paras.gamma[p] << endl;
//                         cout << "gamma size : " << resP->gamma.size() << "," << resP->gamma[p].size() << "," << resP->gamma[p][y_iter].size() << endl;

                        resP->gamma[p][y_iter][x_iter] = modP->paras.gamma[p];

//                         cout << "chi: " << modP->paras.chi[p] << endl;
                        resP->chi[p][y_iter][x_iter] = modP->paras.chi[p];
                }

//                 cout << "store special" << endl;
                if (mode_stats == 1)
                {
//                         cout << "q: " << modP->paras.q[p] << endl;
//                         cout << "q size : " << resP->q.size() << endl;//"," << resP->q[p].size() << "," << resP->q[p][y_iter].size() << endl;
                        resP->q[p][y_iter][x_iter] = modP->paras.q[p];

//                         cout << "alpha_raw" << endl;
                        resP->alpha_raw[p][y_iter][x_iter] = modP->paras.alpha_raw[p];
//                         cout << "alpha" << endl;
                        resP->alpha[p][y_iter][x_iter] = modP->paras.alpha[p];
//                         cout << "sigma_V" << endl;
                        resP->sigma_V[p][y_iter][x_iter] = modP->paras.sigma_V[p];
//                         cout << "I_bal" << endl;
                        resP->I_balance[p][y_iter][x_iter] = modP->paras.I_balance[p];
                }

//                 cout << "done?!" << endl;
                if (mode_stats == 3)
                {
                        resP->q[p][y_iter][x_iter] = modP->paras.q[p];

//                         cout << "KL" << endl;
                        resP->KL_entropy[p][y_iter][x_iter] = modP->paras.KL[p];
//                         cout << "entropy" << endl;
                        resP->entropy[p][y_iter][x_iter] = modP->paras.entropy[p];
                }

                if (mode_stats == 4)
                {
                        resP->KL_entropy[p][y_iter][x_iter] = modP->paras.KL[p];
                        resP->entropy[p][y_iter][x_iter] = modP->paras.entropy[p];
                }

//                 cout << "storing regions" << endl;
                if ((mode_stats == 0) || (mode_stats == 4))
                {
                        resP->regions[p][y_iter][x_iter] = modP->paras.regions[p];
                }
        }
}



void simulation::store_results_approx(simulation *simP, model *modP, results * resP)
{
//         unsigned a = resP->gamma_approx.size() - 1;
//
// //         resP->rate[a].push_back(paras.rate);
        for (int p = 0; p < modP->paras.Npop; p++)
        {
                if (mode_stats == 3)
                        resP->q_approx[p][y_iter][x_iter] = modP->paras.q[p];

//         cout << "gamma: " << modP->paras.gamma[0] << endl;
//         cout << "gamma_size: " << resP->gamma_approx.size() << endl;
                resP->gamma_approx[p][y_iter][x_iter] = modP->paras.gamma[p];
        //         cout << "chi: " << modP->paras.chi[0] << endl;
                resP->chi_approx[p][y_iter][x_iter] = modP->paras.chi[p];
//         resP->regions_approx[a].push_back(modP->paras.regions[0]);
        }
}

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
    //         struct parameters_int *paras = (struct parameters_int *) parameters;

        double rate_max = paras.rate_max;
        double delta = paras.delta;
        double gamma = paras.gamma;

        double rate_ratio = nu/rate_max;

        long double rho_exact = (nu > 0) ? gamma/(rate_max*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(delta)/2)*pow(rate_ratio,gsl_pow_2(gamma)-1)*cosh(gamma*delta*sqrt(-2*log(rate_ratio))) : 0;

//         cout << "rho exact: " << rho_exact << ",\t rho approx: " << rho_approx << endl;

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

        double rate_max = paras.rate_max;
        double delta = paras.delta;
        double gamma = paras.gamma;

        double delta_approx = paras.delta_approx;
        double gamma_approx = paras.gamma_approx;

        double rate_ratio = nu/rate_max;

        long double rho_exact, rho_approx;

        if (nu == 0)
        {
            rho_exact = gamma > 1 ? 0 : 1.0/0.0;
            rho_approx = gamma_approx > 1 ? 0 : 1.0/0.0;
        }
        else
        {
            rho_exact = gamma/(rate_max*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(delta)/2)*pow(rate_ratio,gsl_pow_2(gamma)-1)*cosh(gamma*delta*sqrt(-2*log(rate_ratio)));
            rho_approx = gamma_approx/(rate_max*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(delta_approx)/2)*pow(rate_ratio,gsl_pow_2(gamma_approx)-1)*cosh(gamma_approx*delta_approx*sqrt(-2*log(rate_ratio)));
        }
//         cout << "rho exact: " << rho_exact << ",\t rho approx: " << rho_approx << endl;

        return rho_exact*log2(rho_exact/rho_approx);
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
        gsl_integration_qag(&p_integrate, lower, upper, 1e-7, 1e-7, 10000, 6, ww, &res, &err);		//! integrate G from xth2 to infty
        gsl_integration_workspace_free (ww);

//         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
        return res*(upper-lower);
}

double int_distribution_exact(double nu, void *params)
{
        struct parameters_int paras = *(struct parameters_int *) params;
    //         struct parameters_int *paras = (struct parameters_int *) parameters;

        double rate_max = paras.rate_max;
        double delta = paras.delta;
        double gamma = paras.gamma;

        double rate_ratio = nu/rate_max;
        return (nu > 0) ? gamma/(rate_max*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(delta)/2)*pow(rate_ratio,gsl_pow_2(gamma)-1)*cosh(gamma*delta*sqrt(-2*log(rate_ratio))) : 0;
}




double poisson_distr(int k, double lambda)
{
        return exp(k*log(lambda)-lambda-lgamma(k+1.0));
}


void draw_rates(model *modP, computation *comP, results *resP)
{
//         cout << "start! k: " << comP->k << ", N: " << comP->N << endl;
        resP->rate_inf[comP->k].resize(comP->N);
        // initiate variables
        double rate_inf_tmp;

        // initiate random number generator
        gsl_rng_env_setup();
        const gsl_rng_type *TYPE = gsl_rng_mt19937;
        gsl_rng * rng = gsl_rng_alloc(TYPE);
//         printf("r is a '%s' generator\n", gsl_rng_name(rng));
        gsl_rng_set(rng, 1 + comP->seed_theory[comP->k]);


        // for each neuron get rate_inf
        for (int n=0; n<comP->N; ++n)
        {
                // rejection sampling algorithm
                bool found = false;
                do
                {
                        rate_inf_tmp = gsl_rng_uniform_pos(rng)*modP->paras.rate_max[0];                // generate random rate
//                         cout << "\t nu = " << rate_inf_tmp << " ,\t p(nu) = " << modP->distribution_exact(rate_inf_tmp,0) << endl;
                        // and accept sample according to theoretically predicted probability
                        if (gsl_ran_flat(rng, 0,resP->max_prob) < modP->distribution_exact(rate_inf_tmp,0))
                        {
                                resP->rate_inf[comP->k][n] = rate_inf_tmp;
//                                 resP->N_AP[comP->k][n] = int(rate_inf_tmp*comP->T);
                                found = true;
                        }
                }
                while (found == false);
        }
}

void draw_samples(computation *comP, results *resP)
{
        resP->rate_T[comP->k][comP->j].resize(comP->N,0);
        resP->N_AP[comP->k][comP->j].resize(comP->N,0);

        random_device rd;

//         seed: comP->seed_time[comP->k*comP->draw_finite_time]
//         cout << "random generator test: " << rd << endl;

        mt19937 gen(rd()); // this gives almost deterministic results, for whatever reason...
        for (int n=0; n<comP->N; ++n)
        {
                poisson_distribution<> d(resP->rate_inf[comP->k][n]*comP->T);    // initiating poisson distribution with mean rate_inf*T

                resP->N_AP[comP->k][comP->j][n] = d(gen);
                resP->rate_T[comP->k][comP->j][n] = resP->N_AP[comP->k][comP->j][n]/comP->T;                       // generating a random number of spikes from this
        }
}

// vector<double> get_density_estimate(vector<int> data, computation sim, string kernel)
// {
//         vector<double> p_est(simP->AP_max,0);
//         cout << "getting density estimation..." << endl;
//     //     vector<vector<double> > p_est_single (simP->N);
//         for (int n=0; n<simP->N; ++n)
//         {
//                 if (kernel.compare("poisson") == 0)
//                         for (int i=0; i<simP->AP_max; ++i)
//                                 p_est[i] += poisson_distr(data[n],i)/(simP->N/simP->T);
//                 else
//                         cout << "kernel '" << kernel << "' not yet implemented" << endl;
//         }
//     //     for (int i=0; i<N_max; ++i)
//     //         cout << "prob(i=" << i << "): " << p_est[i] << endl;
//         return p_est;
// }


// void bayesian_estimate(model *modP, computation *comP, results *resP)
// {
//         // wanna calculate: P(lambda|k) = P(k|lambda)*p(lambda)/p(k)
//         // if no prior beliefs exist, p(lambda) = p(k) -> P(lambda|k) = P(k|lambda)
//         // however, if we take p(lambda) to be firing rate distribution from theory, we obtain the following:
//
// //     P_k_lambda  = poisson_distr(k,nu*T);
// //     p_lambda    = modP->distribution_exact(nu,0);
// //     p_k = int(P_k_lambda*p_lambda) d lambda
//
//         resP->p_bayes_est[comP->k][comP->j].resize(resP->steps,0);
//         for (int n=0; n<comP->N; ++n)
//         {
// //                 cout << "APs for neuron " << n << ": " << k << endl;
//                 double p_k = 0;
//
//                 for (int i=1; i < resP->steps; ++i)    // "integration" over lambda    //should also work with integration (faster?)
//                         p_k += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior);
//
//                 resP->p_bayes_est[comP->k][comP->j][0] = 0;
//                 for (int i=1; i < resP->steps; ++i)    // construction of posterior
//                         resP->p_bayes_est[comP->k][comP->j][i] += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior)/p_k / (resP->d_nu*comP->N);
//         }
// }
//
//
// double bayes_est_prior(double nu, model *modP, string prior)
// {
//         if (not prior.compare("mean_field"))
//             return modP->distribution_exact(nu,0);
//         else if (not prior.compare("plain"))
//             return 1;
//         else
//         {
//             cout << "Specify a prior!" << endl;
//             exit(0);
//         }
// }


vector<double> get_cdf(vector<double> p, int steps, double d_nu)     // compute cdf
{
        vector<double> cdf(steps,0);
        for(int i=1; i<=steps; ++i)
                cdf[i] = cdf[i-1] + p[i]*d_nu;

        return cdf;
}

void post_process(computation *comP, model *modP, model *mod_approxP, results *resP)
{
//         resP->KS.resize(comP->draw_from_theory);
//         resP->KL_entropy_single.resize(resP->steps);

//         cout << "gamma: " << modP->paras.gamma[0] << "," << mod_approxP->paras.gamma[0] << endl;
//         cout << "delta: " << modP->paras.delta[0] << "," << mod_approxP->paras.delta[0] << endl;

//         double  KL_tmp = KL_divergence(0,modP->paras.rate_max[0],modP->paras,mod_approxP->paras);
//         cout << "KL = " << KL_tmp << endl;
//         for (int i=1; i<resP->steps; ++i)
//         {
//                 resP->KL_entropy_single[i] = KL_divergence((i-1)*resP->d_nu,i*resP->d_nu,modP->paras,mod_approxP->paras);
//                 cout << "KL[" << i << "] = " << resP->KL_entropy_single[i] << endl;
//         }

        for (int k=0; k<comP->draw_from_theory; ++k)
        {
                resP->KS[k].resize(comP->draw_finite_time,0);
//                 resP->KL_entropy[k].resize(comP->draw_finite_time,0);

                for (int j=0; j<comP->draw_finite_time; ++j)
                {
                        vector<double> cdf_tmp = get_cdf(resP->p_bayes_est[k][j],resP->steps,resP->d_nu);

                        double KS_tmp = 0;
                        for (int i=1; i<resP->steps; ++i)
                        {
                                KS_tmp = abs(cdf_tmp[i] - resP->cdf_theory[i]);
//                                 cout << "KS tmp: " << KS_tmp << ",\tKS old: " << resP->KS[k][j] << endl;
                                resP->KS[k][j] = (KS_tmp > resP->KS[k][j]) ? KS_tmp : resP->KS[k][j];

//                                 cout << "p[" << i << "]: " << resP->p_exact[i] << " vs " << resP->p_bayes_est[k][j][i] << ", log2: " << log2(resP->p_exact[i]/resP->p_bayes_est[k][j][i]) << endl;
//                                 if(resP->p_bayes_est[k][j][i]>pow(10,-3) && resP->p_exact[i]>pow(10,-3))
//                                     resP->KL_entropy[k][j] += resP->p_exact[i] * resP->d_nu * log2(resP->p_exact[i]/resP->p_bayes_est[k][j][i]);
                        }
//                         cout << "KS test: " << resP->KS[k][j] << endl;
                }
        }

}
