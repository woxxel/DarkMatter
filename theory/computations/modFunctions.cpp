#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include "functions.h"

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

    trans_DM.resize(paras.Npop);
    trans_np.resize(paras.Npop);
    trans_imp.resize(paras.Npop);
    trans_inc.resize(paras.Npop);

    trans_DM_found.resize(paras.Npop);
    trans_np_found.resize(paras.Npop);

    in_DM.resize(paras.Npop);
    in_np.resize(paras.Npop);
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
	for (unsigned p = 0; p < paras.Npop; p++)
	{
		// get sigma_V
        // cout << "tau_G: " << paras.tau_G << ", tau_A: " << paras.tau_A << ", tau_N: " << paras.tau_N << ", tau_M: " << paras.tau_M << endl;
        // from excitatory AMPA synapses
        double var_V_A = gsl_pow_2(paras.J_E[p]) * paras.kappa * paras.rate[p] / (paras.tau_A + paras.tau_M) * ( gsl_pow_2(1-paras.n)/2 + (1-paras.n)*paras.n*paras.tau_A / (paras.tau_A + paras.tau_N) );
                // from excitatory NMDA synapses
		double var_V_N = gsl_pow_2(paras.J_E[p]) * paras.kappa * paras.rate[p] / (paras.tau_N + paras.tau_M) * ( paras.n*paras.n/2 + (1-paras.n)*paras.n*paras.tau_N / (paras.tau_A + paras.tau_N) );
                // from inhibitory GABA synapses
		double var_V_G = gsl_pow_2(paras.J_I[p]) * paras.rate[p] * 0.5 / (paras.tau_G + paras.tau_M);

        // from external drive
        // cout << "J_I = " << paras.J_I[p] << ", J_0 = " << paras.J_0 << endl;

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
    for (unsigned p = 0; p < paras.Npop; p++)
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
    for (unsigned p = 0; p < paras.Npop; p++)
    {

        // cout << "nu="<<paras.rate[p]<< ", rate_max=" << paras.rate_max[p] << endl;
        // cout << "alpha="<<paras.alpha[p]<< ", sigma=" << paras.sigma_V[p] << endl;
        // cout << "I: " << I_squared_nu(paras.alpha[p],paras.sigma_V[p],paras.rate[p],paras.rate_max[p]) << endl;
        double I_balance = sqrt( I_squared_nu(paras.alpha[p],paras.sigma_V[p],paras.rate[p],paras.rate_max[p]) );
        paras.I_balance[p] = I_balance;
        paras.delta[p] = I_balance/paras.alpha[p];
        // cout << "got I balance: " << I_balance << endl;
        // cout << "got delta: " << paras.delta[p] << endl;
    }
}

//! all those could have an added boolean, showing whether they were calculated yet, and only evaluate, if not. however, this would require boolean to be updated, whenever some variable is changed
void model::get_gamma()
{
    for (unsigned p = 0; p < paras.Npop; p++)
        paras.gamma[p] = paras.sigma_V[p]/paras.alpha[p];
}

void model::get_chi()
{
    for (unsigned p = 0; p < paras.Npop; p++)
	{
        double nu_peak_log_I = nu_peak_log_full(p);
        paras.chi[p] = -log10(exp(1)) * nu_peak_log_I + log10(paras.rate[p]);
    }
}

double model::nu_peak_log_full(unsigned p)
{
    return log(paras.rate_max[p]) - (gsl_pow_2(paras.gamma[p]*paras.delta[p])-2*(gsl_pow_2(paras.gamma[p])- 1) + paras.gamma[p]*paras.delta[p]*sqrt(gsl_pow_2(paras.gamma[p]*paras.delta[p])-4*(gsl_pow_2(paras.gamma[p])-1))) /(4*gsl_pow_2(gsl_pow_2(paras.gamma[p]) - 1));
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
                for (unsigned p = 0; p < paras.Npop; p++)
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
                for (unsigned p = 0; p < paras.Npop; p++)
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
                for (unsigned p = 0; p < paras.Npop; ++p)
                {
                    double tau_q = 2*(paras.tau_M + paras.tau_G);
//                     cout << "tau_M=" << paras.tau_M << ", tau_G=" << paras.tau_G << endl;
                    double q0 = gsl_pow_2(paras.alpha_0[p])/gsl_pow_2(paras.J_I[p]);
                    double potent = tau_q*(gsl_pow_2(paras.rate[p]) + q0)/(paras.rate[p] + 2*tau_q*(gsl_pow_2(paras.rate[p]) + q0));


                    paras.q[p] = 1./sqrt(1+2*tau_q*(paras.rate[p] + q0/paras.rate[p]))*pow(1+tau_q*(paras.rate[p] + q0/paras.rate[p]),1-potent)*pow(paras.rate_max[p],2*potent)*pow(paras.rate[p],2-2*potent);
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
    bool old_in_inc = in_inc;
    in_inc = false;
    for (unsigned p = 0; p < paras.Npop; p++)
        if (q_border1(p) || q_border2(p))
            in_inc = true || in_inc;

    return in_inc!=old_in_inc;
}

bool model::no_peak(unsigned p)
{
    bool old_in_np = in_np[p];
    in_np[p] = (gsl_pow_2(paras.gamma[p] * paras.delta[p])-4*(gsl_pow_2(paras.gamma[p]) - 1)) < 0;
    return in_np[p]!=old_in_np;
}

bool model::implausible()
{
    bool old_in_imp = in_imp;

    double lower = 0.9;
    double upper = 1;

//         cout << "rate max: " << paras.rate_max[0] << ", gamma: " << paras.gamma[0] << ", delta: " << paras.delta[0] << endl;
    in_imp = false;
    for (unsigned p = 0; p < paras.Npop; p++)
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
        // if (res > 0.1)
        in_imp = res > 0.1 || in_imp;
            // return true;
    }
    return in_imp!=old_in_imp;
}

bool model::DM_border(unsigned p)
{
    bool old_in_DM = in_DM[p];
    in_DM[p] = paras.gamma[p] > 1;
    return in_DM[p]!=old_in_DM;
}

void model::find_transitions(simulation *simP, results *resP)
{

    // cout << "check DM border & no peak" << endl;
    for (unsigned p = 0; p < paras.Npop; p++)
    {
        trans_DM_found[p] = false;
        trans_np_found[p] = false;
        if (DM_border(p) && (simP->vars[0].iter>0)) {
            // if (!trans_DM_found[p]) {
                // cout << " found DM transition for p=" << p << " and iter #" << simP->vars[1].iter << ", @ rate=" << modP->paras.rate[p] << endl;
                // resP->trans_DM[p][simP->vars[1].iter] = paras.rate[p];      // doesnt work in cases, when there is a second DM transition at high nu
            // cout << "DM transition @ " << *simP->vars[0].paraP << endl;
            trans_DM[p] = *simP->vars[0].paraP;
            // trans_DM[p] = paras.rate[p];
            trans_DM_found[p] = true;
            // }
        } else {
            trans_DM_found[p] = false;
        }

        if (no_peak(p) && (simP->vars[0].iter>0)) {
            // if (!trans_np_found[p]) {
                // cout << " found no peak transition for p=" << p << " and iter #" << simP->y_iter << ", @ rate=" << modP->paras.rate[p] << endl;
                // resP->trans_np[p][simP->vars[1].iter] = paras.rate[p];
            trans_np[p] = *simP->vars[0].paraP;
            trans_np_found[p] = true;
            // }
        } else {
        // trans_np_found[p] = false;
        }
    }

    // cout << "check implausible" << endl;
	// starting to find boundaries
    trans_imp_found = false;
	if (implausible() && (simP->vars[0].iter>0)) {
		// if (!trans_imp_found) {
		for (unsigned p = 0; p < paras.Npop; p++)
            // resP->trans_imp[p][simP->vars[1].iter] = paras.rate[p];
            trans_imp[p] = *simP->vars[0].paraP;
		trans_imp_found = true;
		// }
	} else {
        // trans_imp_found = false;
	}

    // cout << "check inconsistent" << endl;
    trans_inc_found = false;
	if (inconsistent() && (simP->vars[0].iter>0)) {
		// if (!trans_inc_found) {
		trans_inc_found = true;
		for (unsigned p = 0; p < paras.Npop; p++)
		{
            // resP->trans_inc[p][simP->vars[1].iter] = paras.rate[p];
			trans_inc[p] = *simP->vars[0].paraP;

			// if (!trans_DM_found[p]) {
			trans_DM[p] = NAN;
			trans_DM_found[p] = true;
			// }

			// if (!trans_np_found[p]) {
				// cout << "forcing np end @ rate = " << modP->paras.rate[p] << endl;
				// cout << "previous value: " << resP->trans_np[p][simP->vars[1].iter]<< endl;
				// resP->trans_np[p][simP->vars[1].iter] = NAN;
			trans_np[p] = NAN;
			trans_np_found[p] = true;
			// }

			// if (!trans_imp_found) {
			trans_imp[p] = NAN;
			trans_imp_found = true;
			// }
		}
		// }
	} else
		trans_inc_found = false;

	if (trans_inc_found && simP->mode_stats==3)
		for (unsigned p = 0; p < paras.Npop; p++)
			resP->KL_entropy[p][simP->vars[1].iter][simP->vars[0].iter] = NAN;

    // cout << "assign region" << endl;
	if ((simP->mode_stats == 0) || (simP->mode_stats == 3))
	{
        int val;
		for (unsigned p = 0; p < paras.Npop; ++p)
		{
			if (in_inc) val = 3;
			else if (in_np[p]) val = 2;
			else if (in_imp) val = 1;
			else val = 0;

			paras.regions[p] = val;
            // regions[simP->n_iter][simP->alpha_0_iter][simP->tau_G_iter][simP->rateWnt_iter][p] = val;
		}
	}

    // cout << "final check" << endl;
	// check if this is the last iteration along x-axis
	if (simP->vars[0].iter==simP->vars[0].steps-1 || in_inc)// && (simP->mode_stats==2))
	{
		for (unsigned p = 0; p < paras.Npop; ++p)
		{
			if (!simP->trans_DM_found[p])
				trans_DM[p] = NAN;
			if (!simP->trans_np_found[p])
				trans_np[p] = NAN;
			if (!simP->trans_imp_found)
				trans_imp[p] = NAN;
			if (!in_inc && !simP->trans_inc_found)
				trans_inc[p] = NAN;
		}
	}
};



void model::integrate_information(info_paras infoParas)
{
    double lower = 0;
    double upper = 0.99;
    // cout << "rate max: " << paras.rate_max[0] << ", gamma: " << paras.gamma[0] << ", delta: " << paras.delta[0] << endl;
    for (unsigned p = 0; p < paras.Npop; p++)
    {
        if (in_inc)
            infoContent[p] = NAN;
        else {
            struct parameters_int Pparas;

            Pparas.rate_max = paras.rate_max[p];
            Pparas.gamma = paras.gamma[p];
            Pparas.delta = paras.delta[p];

            // for now, hardcoded - but change!
            // Pparas.nu0 = infoParas.nu0;
            // Pparas.c = infoParas.c;

            // double dZeta = (infoParas.maxZeta - infoParas.minZeta)/infoParas.nZeta;

            // cout << "zeta: [" << infoParas.minZeta << "," << infoParas.maxZeta << "] , steps: " << infoParas.nZeta << endl;
            // for (unsigned z=0; z<infoParas.nZeta; z++)
            // {
            // Pparas.zeta = infoParas.minZeta + z*dZeta;

            Pparas.I_alpha = paras.I_alpha;
            Pparas.I_beta = paras.I_beta;

            gsl_function p_integrate;
            p_integrate.function = &int_information_distribution; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
            p_integrate.params = &Pparas;
            // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
            gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
            gsl_set_error_handler_off();
            double res, err;
            gsl_integration_qags(&p_integrate, lower*paras.rate_max[p], upper*paras.rate_max[p], 1e-7, 1e-7, 10000, ww, &res, &err);
            gsl_integration_workspace_free (ww);

            // infoContent[p][z] = res;
            infoContent[p] = res;
            // cout << "information content: " << res << endl;
            // }

            // cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
            // if (res > 0.1)

        }

    }
    // return res;
}
