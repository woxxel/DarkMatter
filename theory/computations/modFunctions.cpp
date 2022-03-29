#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include "functions.h"

void PSP::print_PSP()
{
    cout << "\t\t ### Synapse " << s << " properties ###" << endl;
    cout << "\t\t - tau_I: " << tau_I << endl;
    cout << "\t\t - tau_norm: " << tau_norm << endl;
}

void Population::print_population()
{
    cout << "\t ### Population " << p << " properties ###" << endl;
    cout << "\t - I_ext: " << I_ext << endl;
    cout << "\t - rateWnt: " << rateWnt << endl;
    cout << "\t - alpha_0: " << alpha_0 << endl;
    cout << "\t - tau_M: " << tau_M << endl;
    cout << "\t - tau_n: " << tau_n << endl;
    cout << "\t - J0: " << J0 << endl;
    for (unsigned s=0; s<nPSP; s++) {
        psp[s].print_PSP();
    }
}

void Layer::print_layer()
{
    cout << "### Layer " << l << " properties ###" << endl;
    cout << " - eps: " << eps << endl;
    cout << " - eta: " << eta << endl;
    cout << " - J_l: ";
    for (unsigned l=0; l<J_l.size(); l++)
        cout << J_l[l] << ",";
    cout << endl;
    cout << " - kappa: " << kappa << endl;
    for (unsigned p=0; p<nPop; p++) {
        population[p].print_population();
    }
}

// void Model::resize(Simulation *simP)
// {
    // resize all to having p populations

    // paras.pop.resize(paras.Npop);

    // paras.rate.resize(paras.Npop);
    // paras.alpha_0.resize(paras.Npop);
    //
	// paras.J_E.resize(paras.Npop);
	// paras.J_I.resize(paras.Npop);
    // paras.kappa.resize(paras.Npop);


    // parameters to be calculated
    // shouldnt this be in struct "results"?
    // paras.q.resize(paras.Npop);
    //
    // paras.alpha_raw.resize(paras.Npop);
	// paras.alpha.resize(paras.Npop);
	// paras.sigma_V.resize(paras.Npop);
	// paras.rate_max.resize(paras.Npop);
    //
	// paras.gamma.resize(paras.Npop);
	// paras.delta.resize(paras.Npop);
    // paras.I_balance.resize(paras.Npop);
	// paras.chi.resize(paras.Npop);
    // paras.regions.resize(paras.Npop);
    //
    // paras.KL.resize(paras.Npop);
    // paras.entropy.resize(paras.Npop);


    // trans_DM.resize(paras.Npop);
    // trans_np.resize(paras.Npop);
    // trans_imp.resize(paras.Npop);
    // trans_inc.resize(paras.Npop);

    // trans_DM_found.resize(paras.Npop);
    // trans_np_found.resize(paras.Npop);

    // in_DM.resize(paras.Npop);
    // in_np.resize(paras.Npop);
// }

// void Model::add_PSP(int p, double tau_I, double tau_norm, double tau_n)
// {
//     PSP psp = {tau_I, tau_norm, tau_n};
//     paras.pop[p].psp.push_back(psp);
// }

void Model::set_weights()
{
    // implement setter, calling set_weights only when according variables are changed
    for (unsigned l=0; l<L; l++) {
        // if (paras.drive == 2)
        //     paras.J_0 = paras.J * layer[l].population[p].tau_M;

        // set weights for inter-layer coupling
        layer[l].J_l.resize(L);
        for (unsigned ll=0; ll<L; ll++)
            layer[l].J_l[ll] = layer[l].J0_l[ll] * layer[l].population[1].tau_M;

        if (layer[l].nPop == 1)
        {
            layer[l].population[0].J[0] = layer[l].population[0].J0 * layer[l].population[0].tau_M;		   // J_II
            layer[l].population[0].J[1] = layer[l].population[0].J0 * layer[l].population[0].tau_M;		   // J_EI

            layer[l].population[1].J[0] = 0;			       // J_IE
            layer[l].population[1].J[1] = 0;				   // J_EE
        }
        else if (layer[l].nPop == 2)
        {
            // watch out for indexing:
            // inhibitory population: index 0
            // excitatory population: index 1
            layer[l].population[0].J[0] = layer[l].population[0].J0 * sqrt(1 - gsl_pow_2(layer[l].eps)) * layer[l].population[0].tau_M;  // J_II

            layer[l].population[0].J[1] = layer[l].population[0].J0 * sqrt(1 - gsl_pow_2(layer[l].eta * layer[l].eps)) * layer[l].population[0].tau_M;	// J_EI

            layer[l].population[1].J[0] = layer[l].population[1].J0 * layer[l].eps * layer[l].population[1].tau_M;               // J_IE

            layer[l].population[1].J[1] = layer[l].population[1].J0 * layer[l].eta * layer[l].eps * layer[l].population[1].tau_M;   // J_EE

            //! for p=1, the excitatory population receives inhibition, but is decoupled, such that it gives no feedback
            // cout << "eta: " << layer[l].eta << ", eps: " << layer[l].eps << endl;
            // cout << "J_II: " << layer[l].population[0].J[0] << ", J_EI: " << layer[l].population[0].J[1] << endl;
            // cout << "J_IE: " << layer[l].population[1].J[0] << ", J_EE: " << layer[l].population[1].J[1] << endl;
        }
        else
        {
            cout << "There is no algorithm set yet to assign weights for more than 2 populations. Please fix that before you continue!" << endl;
            exit(0);
        }
    }
}

void Model::set_mixture()
{
    for (unsigned l=0; l<L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            if (layer[l].population[p].nPSP==1)
                layer[l].population[p].psp[0].tau_n = 1.;
            else {
                layer[l].population[p].psp[0].tau_n = 1.-layer[l].population[p].tau_n;
                layer[l].population[p].psp[1].tau_n = layer[l].population[p].tau_n;
            }
        }
    }
}

void Model::get_sigma_V()
{
    // get sigma_V

	// iterate over both populations
    double prefactor, J;
    double var_V, var_V_dot, tmp_var_V;
    for (unsigned l=0; l<L; l++) {          // receiving layer
    	for (unsigned p = 0; p < layer[l].nPop; p++) {          // receiving population
            var_V = var_V_dot = 0;  // reset values
            for (unsigned ll=0; ll<L; ll++) {   // projecting layer
            	for (unsigned pp = 0; pp < layer[ll].nPop; pp++) {   // projecting population
                    if (ll!=l && (p==0 || pp==0)) continue;     // inter-layer impact only between excitatory populations

                    if (layer[ll].population[pp].nPSP>2) {
                        cout << "Please specify 1 or 2 types of synapses per population. More are not yet implemented" << endl;
                        throw;
                    }

                    // switch betweeen definitions for inter-layer and intra-layer coupling strengths
                    if (l==ll)  J = layer[ll].population[pp].J[p];
                    else        J = layer[ll].J_l[l];

                    for (unsigned s=0; s<layer[ll].population[pp].nPSP; s++) {

                        // define common prefactor for temporal variance
                        prefactor = gsl_pow_2(J) * layer[ll].population[pp].rateWnt / (layer[ll].population[pp].psp[s].tau_I + layer[l].population[p].tau_M);
                        if (pp==0) prefactor *= layer[ll].kappa;

                        // multiply by factor, specified by populations synaptic transmissions
                        if (layer[ll].population[pp].nPSP==1) tmp_var_V = prefactor / 2.;
                        else if (layer[ll].population[pp].nPSP==2) {

                            tmp_var_V = prefactor * ( gsl_pow_2(layer[ll].population[pp].psp[s].tau_n)/2 + (1-layer[ll].population[pp].psp[s].tau_n)*layer[ll].population[pp].psp[s].tau_n*layer[ll].population[pp].psp[s].tau_I / (layer[ll].population[pp].psp[0].tau_I + layer[ll].population[pp].psp[1].tau_I) );
                        }
                        var_V += tmp_var_V;
                        var_V_dot += tmp_var_V/(layer[ll].population[pp].psp[s].tau_I * layer[l].population[p].tau_M);
                    }
                }
            }
            // if (paras.drive == 2)
            // {
            //     double var_V_0 = sqrt(1/paras.K_0) * paras.J_0 * paras.J_I[p] * paras.rate[p] * 0.5 / (paras.tau_0 + paras.tau_M);
            //     var_V += var_V_0;
            //     var_V_dot += var_V_0 / (paras.tau_0 * paras.tau_M);
            // }
            layer[l].population[p].simulation.sigma_V = sqrt(var_V);

            // and the maximum firing rate response
            layer[l].population[p].simulation.rate_max = sqrt(var_V_dot / var_V) / (2 * M_PI);
            // cout << "(" << l << "," << p << ") ";
            // cout << "eps: " << layer[l].eps << ", mixture: " << layer[l].population[p].psp[0].tau_n << ", var total: " << sqrt(var_V) << ", rate max: " << layer[l].population[p].simulation.rate_max << endl; //", from external sources: " << var_V_0 << endl;
            //                 cout << "sigma population " << p << ": " << paras.sigma_V[p] << endl;
    	}
    }
}

vector< vector<double> > Model::calc_alpha(vector<vector<double>> q)
{
    double J;
    double alpha_sq, alpha_sq_0, tmp_alpha_sq;
    vector< vector<double> > alpha;
    // cout << "calculating alpha" << endl;

    alpha.resize(L);
    for (unsigned l = 0; l < L; l++) {        // receiving layer
        alpha[l].resize(layer[l].nPop);
        for (unsigned p = 0; p < layer[l].nPop; p++) {   // receiving population
            alpha_sq = 0;
            for (unsigned ll = 0; ll < L; ll++) {     // projecting layer
            	for (unsigned pp = 0; pp < layer[ll].nPop; pp++) {   // projecting population
                    if (ll!=l && (p==0 || pp==0)) continue;

                    if (l==ll)  J = layer[ll].population[pp].J[p];
                    else        J = layer[ll].J_l[l];
                    // get alpha
                    tmp_alpha_sq = gsl_pow_2(J) * q[ll][pp];
                    if (pp==0) tmp_alpha_sq *= layer[ll].kappa;

                    alpha_sq += tmp_alpha_sq;
                }
            }

            alpha_sq_0 = 0;
            // if (paraP->drive == 2) {
            //     // quenched variance from afferent, spiking drive (gauss distributed synapse numbers)
            //     alpha_sq_0 = sqrt(1./paraP->K_0) * gsl_pow_2(paraP->J_I[p]) * gsl_pow_2(paraP->rate[p]);
            // }

            // cout << "interresults: "<< alpha_sq << ", " << alpha_sq_0 << ", " << gsl_pow_2(layer[l].population[p].alpha_0) << endl;

            alpha[l][p] = sqrt( alpha_sq + alpha_sq_0 + gsl_pow_2(layer[l].population[p].alpha_0) );
            // cout << "total: " << alpha[l][p] << endl;
        }
    }
    // cout << "calculating alpha done!" << endl;
    return alpha;
}

void Model::get_alpha()
{
    vector<vector<double> > alpha, q;

    q.resize(L);
    for (unsigned l = 0; l < L; l++) {
        q[l].resize(layer[l].nPop);
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            q[l][p] = layer[l].population[p].simulation.q;
        }
    }

    alpha = calc_alpha(q);

    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            // cout << "alpha total: " << gsl_pow_2(alpha[l][p]) << endl;
            layer[l].population[p].simulation.alpha = alpha[l][p];
        }
    }
}

void Model::get_delta()
{
    Population_Simulation *simP;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {

            simP = &layer[l].population[p].simulation;
            // cout << "nu="<<paras.rate[p]<< ", rate_max=" << paras.rate_max[p] << endl;
            // cout << "alpha="<<paras.alpha[p]<< ", sigma=" << paras.sigma_V[p] << endl;
            // cout << "I: " << I_squared_nu(paras.alpha[p],paras.sigma_V[p],paras.rate[p],paras.rate_max[p]) << endl;
            double I_balance = sqrt( I_squared_nu(simP->alpha,simP->sigma_V,layer[l].population[p].rateWnt,simP->rate_max) );
            simP->I_balance = I_balance;
            simP->delta = I_balance/simP->alpha;
            // cout << "got I balance: " << simP->I_balance << endl;
            // cout << "got delta: " << simP->delta << endl;
        }
    }
}

//! all those could have an added boolean, showing whether they were calculated yet, and only evaluate, if not. however, this would require boolean to be updated, whenever some variable is changed
void Model::get_gamma()
{
    Population_Simulation *simP;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            simP = &layer[l].population[p].simulation;

            simP->gamma = simP->sigma_V/simP->alpha;
            // cout << "got gamma: " << simP->gamma << endl;
        }
    }
}

void Model::get_chi()
{
    Population_Simulation *popSimP;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            popSimP = &layer[l].population[p].simulation;
            double nu_peak_log_I = nu_peak_log_full(popSimP);
            popSimP->chi = -log10(exp(1)) * nu_peak_log_I + log10(popSimP->rate);
            // cout << "got chi: " << simP->chi << endl;
        }
    }
}

double Model::nu_peak_log_full(Population_Simulation *popSimP)
{
    // return 0;
    return log(popSimP->rate_max) - (gsl_pow_2(popSimP->gamma*popSimP->delta)-2*(gsl_pow_2(popSimP->gamma)- 1) + popSimP->gamma*popSimP->delta*sqrt(gsl_pow_2(popSimP->gamma*popSimP->delta)-4*(gsl_pow_2(popSimP->gamma)-1))) /(4*gsl_pow_2(gsl_pow_2(popSimP->gamma) - 1));
}

double Model::distribution_exact(double nu, int p)
{
    //! should not be evaluated, if gamma*delta > 200 or so, as cosh -> infty
    double rate_ratio = nu/paras.rate_max[p];
//         cout << "ratio: " << rate_ratio << endl;
//         cout << "log: " << log(rate_ratio) << endl;
//         cout << "cosh: " << cosh(paras.gamma[p]*paras.delta[p]*sqrt(-2*log(rate_ratio))) << endl;
    return paras.gamma[p]/(paras.rate_max[p]*sqrt(-M_PI*log(rate_ratio)))*exp(-gsl_pow_2(paras.delta[p])/2)*pow(rate_ratio,gsl_pow_2(paras.gamma[p])-1)*cosh(paras.gamma[p]*paras.delta[p]*sqrt(-2*log(rate_ratio)));
}

void Model::solve_selfcon(int mode_calc)
{
        // initiate variables
//         set_weights();                             // set parameters
        // layer[0].print_layer();
        // layer[1].print_layer();
        for (unsigned l = 0; l < L; l++) {
            for (unsigned p = 0; p < layer[l].nPop; p++)
                layer[l].population[p].simulation.rate = layer[l].population[p].rateWnt;
        }

        get_sigma_V();

        // cout << "solving the selfconsistency equations with mode " << mode_calc << endl;

        //! solving selfconsistency equations(s)
        if (mode_calc == 0) // exact solution
        {
            // numeric root finding to get exact solution
            // set solver parameter
            gsl_vector *q_guess = gsl_vector_alloc (nPop);
            gsl_multiroot_function F;

//                 int f_tmp = this->selfconsistency_f;
//                 F.f = &f_tmp;
            F.f = &selfconsistency_f;
            F.n = nPop;
            F.params = this;
            unsigned pp = 0;
            for (unsigned l = 0; l < L; l++) {
                for (unsigned p = 0; p < layer[l].nPop; p++) {
                    gsl_vector_set(q_guess, pp, gsl_pow_2(layer[l].population[p].simulation.rate));	         // first guess
                    pp++;
                }
            }
            // initiate the solver (dnewton)
            int status;
            const gsl_multiroot_fsolver_type *Tsolv = gsl_multiroot_fsolver_dnewton;
            gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (Tsolv, nPop);
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
            pp = 0;
            for (unsigned l = 0; l < L; l++) {
                for (unsigned p = 0; p < layer[l].nPop; p++) {
                    layer[l].population[p].simulation.q = gsl_vector_get(s->x,pp);
                    // cout << "rate " << layer[l].population[p].rateWnt << " q(" << pp << ",exact)=" << layer[l].population[p].simulation.q << endl;
                    pp++;
                }
            }
            gsl_multiroot_fsolver_free(s);
        }
        if (mode_calc == 1) {// analytic approximation
            cout << "apply approximation" << endl;
            cout << " should be reviewed - definition of tau_q in 2-pop case not proper anymore!" << endl;
            throw;
            // apply "xyz"-approximation to obtain selfconsistent second moment
    //         for (unsigned l = 0; l < L; l++) {
    //             for (unsigned p = 0; p < layer[l].nPop; p++) {
    //                 double tau_q = 2*(layer[l].population[p].tau_M + paras.tau_G);
    // //                     cout << "tau_M=" << paras.tau_M << ", tau_G=" << paras.tau_G << endl;
    //                 double q0 = gsl_pow_2(layer[l].population[p].alpha_0)/gsl_pow_2(paras.J_I[p]);
    //                 double potent = tau_q*(gsl_pow_2(paras.rate[p]) + q0)/(paras.rate[p] + 2*tau_q*(gsl_pow_2(paras.rate[p]) + q0));
    //
    //
    //                 paras.q[p] = 1./sqrt(1+2*tau_q*(paras.rate[p] + q0/paras.rate[p]))*pow(1+tau_q*(paras.rate[p] + q0/paras.rate[p]),1-potent)*pow(paras.rate_max[p],2*potent)*pow(paras.rate[p],2-2*potent);
    // //                     cout << "q(" << p << ",approx)=" << paras.q[p] << endl;
    //             }
    //         }
        }

        get_alpha();
        get_delta();
        get_gamma();
        get_chi();


        // Population_Simulation *popSimP = &layer[0].population[0].simulation;
        // cout << "rate = " << popSimP->rate << ":\t q=" << popSimP->q << "\t alpha=" << popSimP->alpha << " ,\t sigma=" << popSimP->sigma_V << " ,\t delta=" << popSimP->delta << " ,\t gamma=" << popSimP->gamma << " ,\t chi=" << popSimP->chi << endl;
}

bool Model::q_border1(Population_Simulation *popSimP)
{
    return (1 < popSimP->rate/popSimP->rate_max * sqrt(gsl_pow_2(popSimP->alpha) + gsl_pow_2(popSimP->sigma_V))/popSimP->sigma_V);
}

bool Model::q_border2(Population_Simulation *popSimP)
{
    return (1 < popSimP->q/gsl_pow_2(popSimP->rate_max) * sqrt(2 * gsl_pow_2(popSimP->alpha) + gsl_pow_2(popSimP->sigma_V))/popSimP->sigma_V);
}

bool Model::DM_border(Population_Simulation *popSimP)
{
    bool old_in_DM = popSimP->in_DM;
    popSimP->in_DM = popSimP->gamma > 1;
    return popSimP->in_DM!=old_in_DM;
}

bool Model::no_peak(Population_Simulation *popSimP)
{
    bool old_in_np = popSimP->in_np;
    popSimP->in_np = (gsl_pow_2(popSimP->gamma * popSimP->delta)-4*(gsl_pow_2(popSimP->gamma) - 1)) < 0;
    return popSimP->in_np!=old_in_np;
}

bool Model::implausible(Model_Simulation *mSimP)
{
    bool old_in_imp = mSimP->in_imp;

    double lower = 0.9;
    double upper = 1;

//         cout << "rate max: " << paras.rate_max[0] << ", gamma: " << paras.gamma[0] << ", delta: " << paras.delta[0] << endl;
    mSimP->in_imp = false;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {

            struct parameters_int Pparas;

            Pparas.rate_max = layer[l].population[p].simulation.rate_max;
            Pparas.gamma = layer[l].population[p].simulation.gamma;
            Pparas.delta = layer[l].population[p].simulation.delta;

            gsl_function p_integrate;
            p_integrate.function = &int_distribution_exact; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
            p_integrate.params = &Pparas;
            // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
            gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
            gsl_set_error_handler_off();
            double res, err;
    //         try
    //         {
            gsl_integration_qags(&p_integrate, lower*layer[l].population[p].simulation.rate_max, upper*layer[l].population[p].simulation.rate_max, 1e-7, 1e-7, 10000, ww, &res, &err);		//! integrate G from xth2 to infty
    //         } catch (const std::exception& e){
    //             cout << "whatup?" << endl;
    //         } catch (...)
    //         {
    //             cout << "error while integrating..." << endl;
    //         }
            gsl_integration_workspace_free (ww);

    //         cout << "result of integrating distribution over [" << lower << "," << upper << "]: " << res << endl;
            // if (res > 0.1)
            mSimP->in_imp = res > 0.1 || mSimP->in_imp;
            // return true;
        }
    }
    return mSimP->in_imp!=old_in_imp;
}

bool Model::inconsistent(Model_Simulation *mSimP)
{
    bool old_in_inc = mSimP->in_inc;
    mSimP->in_inc = false;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            if (q_border1(&layer[l].population[p].simulation) || q_border2(&layer[l].population[p].simulation))
                mSimP->in_inc = true || mSimP->in_inc;
        }
    }

    return mSimP->in_inc!=old_in_inc;
}

void Model::find_transitions(Simulation *simP)
{
    // cout << "check DM border & no peak" << endl;
    Population_Simulation *popSimP;

    Model_Simulation *mSimP = &simulation;

    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {
            popSimP = &layer[l].population[p].simulation;

            popSimP->trans_DM_found = false;
            popSimP->trans_np_found = false;
            if (DM_border(popSimP) && (simP->vars[0].iter>0)) {
                // if (!trans_DM_found[p]) {
                    // resP->trans_DM[p][simP->vars[1].iter] = paras.rate[p];      // doesnt work in cases, when there is a second DM transition at high nu
                // cout << "DM transition @ " << *simP->vars[0].paraP << endl;
                popSimP->trans_DM = *simP->vars[0].paraP[0];
                // trans_DM[p] = paras.rate[p];
                popSimP->trans_DM_found = true;
                // }
            } else {
                popSimP->trans_DM_found = false;
            }

            if (no_peak(popSimP) && (simP->vars[0].iter>0)) {
                // if (!trans_np_found[p]) {
                    // cout << " found no peak transition for p=" << p << " and iter #" << simP->y_iter << ", @ rate=" << modP->paras.rate[p] << endl;
                    // resP->trans_np[p][simP->vars[1].iter] = paras.rate[p];
                popSimP->trans_np = *simP->vars[0].paraP[0];
                popSimP->trans_np_found = true;
                // }
            }
        }
    }

    // cout << "check implausible" << endl;
	// starting to find boundaries
    // cout << "checking transition for iter #" << simP->vars[1].iter << "," << simP->vars[0].iter << ", @ " << simP->vars[0].name << "=" << *simP->vars[0].paraP << endl;
    mSimP->trans_imp_found = false;
	if (implausible(mSimP) && (simP->vars[0].iter>0)) {
		// if (!trans_imp_found) {
        mSimP->trans_imp = *simP->vars[0].paraP[0];
        // for (unsigned l = 0; l < L; l++) {
            // for (unsigned p = 0; p < layer[l].nPop; p++)
                // resP->trans_imp[p][simP->vars[1].iter] = paras.rate[p];
        // }
		mSimP->trans_imp_found = true;
	}

    // cout << "check inconsistent" << endl;
    mSimP->trans_inc_found = false;
	if (inconsistent(mSimP) && (simP->vars[0].iter>0)) {
		// if (!trans_inc_found) {
		mSimP->trans_inc_found = true;
        mSimP->trans_inc = *simP->vars[0].paraP[0];

        mSimP->trans_imp = NAN;
        mSimP->trans_imp_found = true;

        for (unsigned l = 0; l < L; l++) {
            for (unsigned p = 0; p < layer[l].nPop; p++) {
                popSimP = &layer[l].population[p].simulation;
                // popResP = &layer[l].population[p].results;
                // resP->trans_inc[p][simP->vars[1].iter] = paras.rate[p];

    			// if (!trans_DM_found[p]) {
    			popSimP->trans_DM = NAN;
    			popSimP->trans_DM_found = true;
    			// }

    			// if (!trans_np_found[p]) {
    				// cout << "forcing np end @ rate = " << modP->paras.rate[p] << endl;
    				// cout << "previous value: " << resP->trans_np[p][simP->vars[1].iter]<< endl;
    				// resP->trans_np[p][simP->vars[1].iter] = NAN;
    			popSimP->trans_np = NAN;
    			popSimP->trans_np_found = true;
    			// }

    			// if (!trans_imp_found) {
    			// }
    		}
        }
		// }
	}

	// if (trans_inc_found && simP->mode_stats==3)
	// 	for (unsigned p = 0; p < paras.Npop; p++)
	// 		resP->KL_entropy[p][simP->vars[1].iter][simP->vars[0].iter] = NAN;

    // cout << "assign region" << endl;
	if ((simP->mode_stats == 0) || (simP->mode_stats == 3))
	{
        int val;
        for (unsigned l = 0; l < L; l++) {
            for (unsigned p = 0; p < layer[l].nPop; p++) {
    			if (mSimP->in_inc) val = 3;
    			else if (layer[l].population[p].simulation.in_np) val = 2;
    			else if (mSimP->in_imp) val = 1;
    			else val = 0;

    			layer[l].population[p].simulation.regions = val;
                // regions[simP->n_iter][simP->alpha_0_iter][simP->tau_G_iter][simP->rateWnt_iter][p] = val;
    		}
        }
	}

    // cout << "final check" << endl;
	// check if this is the last iteration along x-axis
	if (simP->vars[0].iter==simP->vars[0].steps-1 || mSimP->in_inc)// && (simP->mode_stats==2))
	{
        for (unsigned l = 0; l < L; l++) {
            for (unsigned p = 0; p < layer[l].nPop; p++) {
                popSimP = &layer[l].population[p].simulation;

    			if (!popSimP->trans_DM_found)
    				popSimP->trans_DM = NAN;
    			if (!popSimP->trans_np_found)
    				popSimP->trans_np = NAN;
    			if (!mSimP->trans_imp_found)
    				mSimP->trans_imp = NAN;
    			if (!mSimP->in_inc && !mSimP->trans_inc_found)
    				mSimP->trans_inc = NAN;
    		}
        }
	}
};



void Model::integrate_information()
{
    Population_Simulation *popSimP;

    double lower = 0;
    double upper = 0.9;
    // cout << "rate max: " << paras.rate_max[0] << ", gamma: " << paras.gamma[0] << ", delta: " << paras.delta[0] << endl;
    for (unsigned l = 0; l < L; l++) {
        for (unsigned p = 0; p < layer[l].nPop; p++) {

            popSimP = &layer[l].population[p].simulation;

            struct parameters_int Pparas;

            Pparas.rate_max = popSimP->rate_max;
            Pparas.gamma = popSimP->gamma;
            Pparas.delta = popSimP->delta;

            // for now, hardcoded - but change!
            // Pparas.nu0 = infoParas.nu0;
            // Pparas.c = infoParas.c;

            // double dZeta = (infoParas.maxZeta - infoParas.minZeta)/infoParas.nZeta;

            // cout << "zeta: [" << infoParas.minZeta << "," << infoParas.maxZeta << "] , steps: " << infoParas.nZeta << endl;
            // for (unsigned z=0; z<infoParas.nZeta; z++)
            // {
            // Pparas.zeta = infoParas.minZeta + z*dZeta;

            Pparas.I_alpha = I_alpha;
            Pparas.I_beta = I_beta;

            gsl_function p_integrate;
            p_integrate.function = &int_information_distribution; //this should be changed to general function (can't be, as this function here needs to have special variables -> integration)
            p_integrate.params = &Pparas;
            // integration method, where function is divided into subintervals. at each iteration, intervall with highest error is split up -> good approximation of function
            gsl_integration_workspace *ww = gsl_integration_workspace_alloc(10000);	//! workspace for storage of integration data (up to 1000 intervalls)
            gsl_set_error_handler_off();
            double res, err;
            gsl_integration_qags(&p_integrate, lower*popSimP->rate_max, upper*popSimP->rate_max, 1e-7, 1e-7, 10000, ww, &res, &err);
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
