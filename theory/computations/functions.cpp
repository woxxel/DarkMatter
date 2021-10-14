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
    simP->trans_DM_found.resize(modP->paras.Npop,false);
    simP->trans_np_found.resize(modP->paras.Npop,false);
    simP->trans_imp_found = false;
    simP->trans_inc_found = false;

    // now: size of border vectors
    resP->trans_DM.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_np.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_inc.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));
    resP->trans_imp.resize(modP->paras.Npop,vector<double>(simP->vars[0].steps));

    resP->rate.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->q.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->gamma.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    resP->chi.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

    // if ((simP->mode_stats == 0) || (simP->mode_stats == 4))
    resP->regions.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

    if (simP->mode_stats == 1)
    {
        // cout << "start q size : " << resP->q.size() << endl;
        resP->alpha_raw.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->alpha.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->sigma_V.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->I_balance.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
    }

    if ((simP->mode_stats == 2) || (simP->mode_stats == 4))
    {
        resP->q_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->steps,vector<double>(simP->steps)));
        resP->gamma_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->chi_approx.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

        resP->KL_entropy.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));
        resP->entropy.resize(modP->paras.Npop,vector<vector<double> >(simP->vars[1].steps,vector<double>(simP->vars[0].steps)));

        resP->p_range.resize(simP->steps);
        resP->p_exact.resize(simP->steps);
        resP->cdf_theory.resize(simP->steps);
    }
}

void compare_approx(model *modP, model *mod_approxP, simulation *simP, results *resP)
{
    for (int p = 0; p < modP->paras.Npop; p++)
    {
        long double cosh_tmp = cosh(modP->paras.gamma[p]*modP->paras.delta[p]*sqrt(-2*log(pow(10,-6)/modP->paras.rate_max[p]))); // check if pdf can be calculated
        if (!isinf(cosh_tmp) && !isnan(cosh_tmp))
        {
            // cout << resP->KL_entropy.size() << endl;
            modP->paras.KL[p] = KL_divergence(p,0,modP->paras.rate_max[p],modP->paras,mod_approxP->paras);
            modP->paras.entropy[p] = shannon_entropy(p,0,modP->paras.rate_max[p],modP->paras);
            // if (simP->steps == 1)
            // {
            // cout << "Calculate pdf and cdf from selfconsistent solution..." << endl;

            resP->d_nu = modP->paras.rate_max[p]/simP->steps;

            resP->p_exact[0] = 0;
            resP->cdf_theory[0] = 0;
//                                                                         cout << "write p for " << simP->steps << " steps with d_nu = " << resP->d_nu << "." << endl;
            for (unsigned i=1; i<simP->steps; ++i)
            {
                resP->p_range[i] = i*resP->d_nu;    // should be written in some seperate function
                resP->p_exact[i] = modP->distribution_exact(resP->p_range[i],0);
//                                                                                     resP->max_prob = max(resP->p_exact[i],resP->max_prob);
                resP->cdf_theory[i] = resP->cdf_theory[i-1] + pdf2hist((i-1)*resP->d_nu,i*resP->d_nu,modP->paras);

                // obtain maximum value of probability density (with bins of width=rate_max/N_max)
            }
            // }
        }
        else
        {
            modP->paras.KL[p] = NAN;
            modP->paras.entropy[p] = NAN;
        }
    }
}

void find_transitions(model *modP, simulation *simP, results *resP)
{
    for (int p = 0; p < modP->paras.Npop; p++)
    {

        if (modP->paras.gamma[p] > 1) {
            if (!simP->trans_DM_found[p]) {
                // cout << " found DM transition for p=" << p << " and iter #" << simP->vars[1].iter << ", @ rate=" << modP->paras.rate[p] << endl;
                resP->trans_DM[p][simP->vars[1].iter] = modP->paras.rate[p];      // doesnt work in cases, when there is a second DM transition at high nu
                simP->trans_DM_found[p] = true;
            }
        } else {
            simP->trans_DM_found[p] = false;
        }

        if (modP->no_peak(p)) {
            if (!simP->trans_np_found[p]) {
                // cout << " found no peak transition for p=" << p << " and iter #" << simP->y_iter << ", @ rate=" << modP->paras.rate[p] << endl;
                resP->trans_np[p][simP->vars[1].iter] = modP->paras.rate[p];
                simP->trans_np_found[p] = true;
            }
        } else {
            // simP->trans_np_found[p] = false;
        }
    }

	// starting to find boundaries
	if (modP->implausible()) {
		if (!simP->trans_imp_found) {
			for (int p = 0; p < modP->paras.Npop; p++)
				resP->trans_imp[p][simP->vars[1].iter] = modP->paras.rate[p];
			simP->trans_imp_found = true;
		}
	} else {
		// simP->trans_imp_found = false;
	}

	if (modP->inconsistent()) {
		if (!simP->trans_inc_found) {
			simP->trans_inc_found = true;
			for (int p = 0; p < modP->paras.Npop; p++)
			{
				resP->trans_inc[p][simP->vars[1].iter] = modP->paras.rate[p];

				if (!simP->trans_DM_found[p]) {
					resP->trans_DM[p][simP->vars[1].iter] = NAN;
					simP->trans_DM_found[p] = true;
				}

				if (!simP->trans_np_found[p]) {
					// cout << "forcing np end @ rate = " << modP->paras.rate[p] << endl;
					// cout << "previous value: " << resP->trans_np[p][simP->vars[1].iter]<< endl;
					resP->trans_np[p][simP->vars[1].iter] = NAN;
					simP->trans_np_found[p] = true;
				}

				if (!simP->trans_imp_found) {
					resP->trans_imp[p][simP->vars[1].iter] = NAN;
					simP->trans_imp_found = true;
				}
			}
		}
	} else {
		simP->trans_inc_found = false;
	}
		//                                                         }
	if (simP->trans_inc_found && simP->mode_stats==3)
		for (int p = 0; p < modP->paras.Npop; p++)
			resP->KL_entropy[p][simP->vars[1].iter][simP->vars[0].iter] = NAN;

//                                     cout << "stuff done" << endl;
	if ((simP->mode_stats == 0) || (simP->mode_stats == 4))
	{
        int val;
		for (int p = 0; p < modP->paras.Npop; ++p)
		{
			if (simP->trans_inc_found) val = 3;
			else if (simP->trans_np_found[p]) val = 2;
			else if (simP->trans_imp_found) val = 1;
			else val = 0;

			modP->paras.regions[p] = val;
//                                                     regions[simP->n_iter][simP->alpha_0_iter][simP->tau_G_iter][simP->rateWnt_iter][p] = val;
		}
	}

	// check if this is the last iteration along x-axis
	if ((simP->vars[0].iter==simP->vars[0].steps-1) && (simP->mode_stats==2))
	{
		for (int p = 0; p < modP->paras.Npop; ++p)
		{
			if (!simP->trans_DM_found[p])
				resP->trans_DM[p][simP->vars[1].iter] = NAN;
			if (!simP->trans_np_found[p])
				resP->trans_np[p][simP->vars[1].iter] = NAN;
			if (!simP->trans_imp_found)
				resP->trans_imp[p][simP->vars[1].iter] = NAN;
			if (!simP->trans_inc_found)
				resP->trans_inc[p][simP->vars[1].iter] = NAN;
		}
	}
};

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
    // implement setter, calling set_weights only when according variables are changed
    if (paras.drive == 2)
        paras.J_0 = paras.J * paras.tau_M;

    if (paras.Npop == 1)
    {
        paras.J_I[0] = paras.J * paras.tau_M;		   // J_II
        paras.J_I[1] = paras.J * paras.tau_M;		   // J_EI

        paras.J_E[0] = 0;			       // J_IE
        paras.J_E[1] = 0;				   // J_EE
        // cout << "new weights: " << paras.J_I[0] << endl;
        // cout << "J: " << paras.J << endl;
    }
    else if (paras.Npop == 2)
    {
        // watch out for indexing:
        // inhibitory population: index 0
        // excitatory population: index 1
        paras.J_I[0] = paras.J * sqrt(1 - gsl_pow_2(paras.eps)) * paras.tau_M;				// J_II
        paras.J_I[1] = paras.J * sqrt(1 - gsl_pow_2(paras.eta * paras.eps)) * paras.tau_M;	// J_EI

        paras.J_E[0] = paras.J * paras.eps * paras.tau_M;			                            // J_IE
        paras.J_E[1] = paras.J * paras.eta * paras.eps * paras.tau_M;					        // J_EE

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



// double simulation::y_val()
// {
//     if (tau_GSz > 1)
//         return tau_G[tau_G_iter];
//     else if (alpha_0Sz > 1)
//         return alpha_0[alpha_0_iter];
//     else if (nSz > 1)
//         return n[n_iter];
//     else
//         return 0;
// }


void simulation::store_results(simulation *simP, model * modP, model * mod_approxP, results * resP)
{
//         unsigned a = resP->rate.size() - 1;
//         cout << "size = " << a << endl;
//         unsigned size = resP->rate[a].size()+1;
//         cout << "rate=" << paras.rate << " ,\t q=" << paras.q[0] << " ,\t alpha=" << paras.alpha[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
        //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent

    for (int p = 0; p < modP->paras.Npop; p++)
    {
        resP->rate[p][vars[1].iter][vars[0].iter] = modP->paras.rate[p];
        resP->q[p][vars[1].iter][vars[0].iter] = modP->paras.q[p];

        resP->gamma[p][vars[1].iter][vars[0].iter] = modP->paras.gamma[p];
        resP->chi[p][vars[1].iter][vars[0].iter] = modP->paras.chi[p];

        resP->regions[p][vars[1].iter][vars[0].iter] = modP->paras.regions[p];

        if (mode_stats == 1)
        {
            resP->alpha_raw[p][vars[1].iter][vars[0].iter] = modP->paras.alpha_raw[p];
            resP->alpha[p][vars[1].iter][vars[0].iter] = modP->paras.alpha[p];
            resP->sigma_V[p][vars[1].iter][vars[0].iter] = modP->paras.sigma_V[p];
            resP->I_balance[p][vars[1].iter][vars[0].iter] = modP->paras.I_balance[p];
        }

        if (mode_stats == 2)
        {
            resP->q_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.q[p];
            resP->gamma_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.gamma[p];
            resP->chi_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.chi[p];

            resP->KL_entropy[p][vars[1].iter][vars[0].iter] = modP->paras.KL[p];
            resP->entropy[p][vars[1].iter][vars[0].iter] = modP->paras.entropy[p];
        }

    }
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
