#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

// #include "spdlog/spdlog.h"
// spdlog::set_level(spdlog::level::debug);

#include "./theory/structures.h"
#include "./theory/functions.cpp"
#include "./io/readwrite.cpp"

// using namespace std;

int main(int argc, char** argv)
{
	if ((argc < 4) or (argc > 4))
	{
		cout << "Please specify (only) the input and output file name!" << endl;
		return 1;
	}

//      parameters:
//          drive:
//          0 : all populations driven to have same firing rates
//          1 : only excitatory one driven
//          2 : recurrent, inhibitory population is driven by afferent spikes from excitatory population
	string fileModel = argv[1];
    string fileSim = argv[2];
	string fileOut = argv[3];

    simulation sim, *simP = &sim;
    model mod, *modP = &mod;
    results res, *resP = &res;

    read_model(fileModel,modP);

    read_simulation(fileSim,simP);

    initiate_results(mod,simP,resP);
	// cout << "sizes (after sim) - n: " << simP->nSz << ", alpha_0: " << simP->alpha_0Sz << ", tau_G: " << simP->tau_GSz << ", rate: " << simP->rateWntSz << endl;

    model mod_approx, *mod_approxP = &mod_approx;
    mod_approx = mod;

    // cout << "calculate solutions for mode stats = " << sim.mode_stats << endl;
    int val;

//         cout << "timeconstants: tau_M = " << mod.paras.tau_M << ", tau_G = " << mod.paras.tau_G << ", tau_A = " << mod.paras.tau_A << ", tau_N = " << mod.paras.tau_N << endl;

    // cout << "steps:" << sim.steps << endl;
    if (sim.steps == 1) sim.initiate_y_axis(mod);

    for (sim.alpha_0_iter = 0; sim.alpha_0_iter < sim.alpha_0Sz; sim.alpha_0_iter++)
    {
        // spdlog::debug("in alpha");
        for (int p = 0; p < mod.paras.Npop; ++p)
		{
			mod.paras.alpha_0[p] = sim.alpha_0[sim.alpha_0_iter];
		}
		// cout << "still here " << endl;

        if (res.axes[0] == 2)
            sim.initiate_y_axis(mod);
		// cout << "and yet, still here " << endl;
		// cout << "n iteration: " << sim.n_iter << "/" << sim.nSz << endl;

        for (sim.n_iter = 0; sim.n_iter < sim.nSz; sim.n_iter++)
        {
            // spdlog::debug("in n");
            mod.paras.n = sim.n[sim.n_iter];

            if (res.axes[1] == 2)
                sim.initiate_y_axis(mod);


            for (sim.eps_iter = 0; sim.eps_iter < sim.epsSz; ++sim.eps_iter)
            {
                // spdlog::debug("in eps");
                // cout << "iter: " << sim.eps_iter << ", value: " << sim.eps[sim.eps_iter] << endl;
                mod.paras.eps = sim.eps[sim.eps_iter];
                if (res.axes[2] == 2)
                        sim.initiate_y_axis(mod);


                for (sim.eta_iter = 0; sim.eta_iter < sim.etaSz; ++sim.eta_iter)
                {
                    // spdlog::debug("in eta");
                    mod.paras.eta = sim.eta[sim.eta_iter];

                    mod.set_weights();  // weights are only influenced by eps, eta and tau_M
                    if (res.axes[3] == 2)
                            sim.initiate_y_axis(mod);


                    for (sim.tau_G_iter = 0; sim.tau_G_iter < sim.tau_GSz; sim.tau_G_iter++)
                    {
                        // spdlog::debug("in tau_G");
                        mod.paras.tau_G = sim.tau_G[sim.tau_G_iter];

                        // cout << "timeconstants: tau_M = " << mod.paras.tau_M << ", tau_G = " << mod.paras.tau_G << ", tau_A = " << mod.paras.tau_A << ", tau_N = " << mod.paras.tau_N << endl;

                        if (mod.paras.tau_G==0)
                            continue;


                        if (res.axes[4] == 2)
                            sim.initiate_y_axis(mod);

                        for (sim.rateWnt_iter = 0; sim.rateWnt_iter < sim.rateWntSz; sim.rateWnt_iter++)
                        {
                            sim.x_iter++;
//                                     if (sim.mode_stats == 1)// && (sim.max_ax[1] == 1))
//                                             sim.y_iter = 0;

                            for (int p = 0; p < mod.paras.Npop; p++)
                            {
                                mod.paras.rate[p] = sim.rateWnt[sim.rateWnt_iter];

                                if (mod.paras.rate[p]==0)
                                    continue;
                            }

//                                     mod.paras.tau_G = mod.paras.rate[0]*(0.030+mod.paras.tau_M) - mod.paras.tau_M;
//                                     cout << "rate: " << mod.paras.rate[0] << ", tau_M: " << mod.paras.tau_M << ", tau_I: " << mod.paras.tau_G << endl;


                            if (mod.paras.drive == 1) // inhibitory rate determined by network parameters
                            {
                                mod.paras.rate[0] = sim.rateWnt[sim.rateWnt_iter]*mod.paras.kappa*mod.paras.J_E[0]/mod.paras.J_I[0];
//                                                             cout << " rate inhib: " << mod.paras.rate[0] << "\t rate exc: " << mod.paras.rate[1] << endl;
                            }
//                                                         cout << "test" << endl;
//                                         cout << "nu: " << sim.rateWnt[sim.rateWnt_iter] << " ,\t tau_G: " << sim.tau_G[sim.tau_G_iter] << " ,\t alpha_0: " << sim.alpha_0[sim.alpha_0_iter] << endl;

//                                                         if (!sim.trans_inc_found)
//                                                         {
//                                                             cout << "Now solve..." << endl;
                                mod.solve_selfcon(sim.mode_calc);
//                                                             cout << "done!" << endl;

                                if ((sim.mode_stats == 3) || (sim.mode_stats == 4))
                                {
                                    mod_approx = mod;
                                    mod_approx.solve_selfcon(1);
                                }


                                for (int p = 0; p < mod.paras.Npop; p++)
                                {
									// cout << "gamma : " << mod.paras.gamma[p] << endl;

                                    if ((sim.mode_stats == 3) || (sim.mode_stats == 4))
                                    {
                                        long double cosh_tmp = cosh(mod.paras.gamma[p]*mod.paras.delta[p]*sqrt(-2*log(pow(10,-6)/mod.paras.rate_max[p]))); // check if pdf can be calculated
                                        if (!isinf(cosh_tmp) && !isnan(cosh_tmp))
                                        {
                                            cout << res.KL_entropy.size() << endl;
                                            mod.paras.KL[p] = KL_divergence(p,0,modP->paras.rate_max[p],modP->paras,mod_approxP->paras);
                                            mod.paras.entropy[p] = shannon_entropy(p,0,modP->paras.rate_max[p],modP->paras);
                                            if (sim.steps == 1)
                                            {
                                                cout << "Calculate pdf and cdf from selfconsistent solution..." << endl;

                                                res.d_nu = modP->paras.rate_max[p]/res.steps;

                                                res.p_exact[0] = 0;
                                                res.cdf_theory[0] = 0;
//                                                                         cout << "write p for " << res.steps << " steps with d_nu = " << res.d_nu << "." << endl;
                                                for (int i=1; i<res.steps; ++i)
                                                {
                                                    res.p_range[i] = i*res.d_nu;    // should be written in some seperate function
                                                    res.p_exact[i] = mod.distribution_exact(res.p_range[i],0);
//                                                                                     res.max_prob = max(res.p_exact[i],res.max_prob);
                                                    res.cdf_theory[i] = res.cdf_theory[i-1] + pdf2hist((i-1)*res.d_nu,i*res.d_nu,mod.paras);

                                                    // obtain maximum value of probability density (with bins of width=rate_max/N_max)
                                                }
                                            }
                                        }
                                        else
                                        {
                                                mod.paras.KL[p] = NAN;
                                                mod.paras.entropy[p] = NAN;
                                        }
                                    }

                                    if ((sim.mode_stats == 0) || (sim.mode_stats == 2) || (sim.mode_stats == 4))
                                    {
//                                                                             if (!sim.trans_DM_found[p])
                                        if (mod.paras.gamma[p] > 1)
                                        {
                                            res.trans_DM[p][sim.y_iter] = mod.paras.rate[p];      // doesnt work in cases, when there is a second DM transition at high nu
                                            sim.trans_DM_found[p] = true;
                                        }
                                        else
                                            sim.trans_DM_found[p] = false;

//                                                                             if (!sim.trans_np_found[p])
                                        if (mod.no_peak(p))
                                        {
                                                res.trans_np[p][sim.y_iter] = mod.paras.rate[p];
                                            sim.trans_np_found[p] = true;
                                        }
                                        else
                                            sim.trans_np_found[p] = false;
                                    }
                                }

                                //                                                             if (!sim.trans_imp_found)
                                if (mod.implausible())
                                {
                                    for (int p = 0; p < mod.paras.Npop; p++)
                                        res.trans_imp[p][sim.y_iter] = mod.paras.rate[p];
                                    sim.trans_imp_found = true;
                                }
                                else
                                    sim.trans_imp_found = false;

                                if (mod.inconsistent())
                                {
                                    sim.trans_inc_found = true;
                                    if ((sim.mode_stats==0) || (sim.mode_stats == 2) || (sim.mode_stats == 4))
                                    {
                                        for (int p = 0; p < mod.paras.Npop; p++)
                                        {
                                            res.trans_inc[p][sim.y_iter] = mod.paras.rate[p];

                                            if (!sim.trans_DM_found[p])
                                                    res.trans_DM[p][sim.y_iter] = NAN;
                                            sim.trans_DM_found[p] = true;

                                            if (!sim.trans_np_found[p])
                                                    res.trans_np[p][sim.y_iter] = NAN;
                                            sim.trans_np_found[p] = true;

                                            if (!sim.trans_imp_found)
                                                    res.trans_imp[p][sim.y_iter] = NAN;
                                            sim.trans_imp_found = true;
//                                                     break;
                                        }
                                    }
                                }
                                else
                                    sim.trans_inc_found = false;
                                //                                                         }
                            if (sim.trans_inc_found && sim.mode_stats==3)
                                for (int p = 0; p < mod.paras.Npop; p++)
                                    res.KL_entropy[p][sim.y_iter][sim.x_iter] = NAN;

//                                     cout << "stuff done" << endl;

                            if ((sim.mode_stats == 0) || (sim.mode_stats == 4))
                            {
                                for (int p = 0; p < mod.paras.Npop; ++p)
                                {
                                    if (sim.trans_inc_found) val = 3;
                                    else if (sim.trans_np_found[p]) val = 2;
                                    else if (sim.trans_imp_found) val = 1;
                                    else val = 0;

                                    mod.paras.regions[p] = val;
//                                                     regions[sim.n_iter][sim.alpha_0_iter][sim.tau_G_iter][sim.rateWnt_iter][p] = val;
                                }
                            }

                            // cout << "storing ..." << endl;
                            sim.store_results(sim, modP,resP);

                            // cout << "storing #2 ..." << endl;
                            if ((sim.mode_stats==3) || (sim.mode_stats==4))
                                sim.store_results_approx(sim, mod_approxP,resP);
//                                                 cout << "KL = " << res.KL_entropy[sim.alpha_0_iter][sim.rateWnt_iter] << endl;
//
//                                                         cout << " done." << endl;

                            // check if this is the last iteration along x-axis
                            if ((sim.rateWnt_iter==sim.steps-1) && (sim.mode_stats==2))
                            {
                                for (int p = 0; p < mod.paras.Npop; ++p)
                                {
                                    if (!sim.trans_DM_found[p])
                                        res.trans_DM[p][sim.y_iter] = NAN;
                                    if (!sim.trans_np_found[p])
                                        res.trans_np[p][sim.y_iter] = NAN;
                                    if (!sim.trans_imp_found)
                                        res.trans_imp[p][sim.y_iter] = NAN;
                                    if (!sim.trans_inc_found)
                                        res.trans_inc[p][sim.y_iter] = NAN;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

// 	for (unsigned i=0;i<sim.steps;++i)
	if (sim.steps == 1)
	{
		write_theory(fileOut,resP);
	} else {
        if ((sim.mode_stats == 0) || (sim.mode_stats == 4))
            write_sharks(fileOut,sim,mod,resP);//,gammaP,chiP,regionsP);
        else
            write_stats(fileOut,mod,sim,resP);
    }
    // cout << "end" << endl;
	return 0;
}
