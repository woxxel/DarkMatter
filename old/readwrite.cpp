#include <iostream>
// #include <netcdfcpp.h>
#include <netcdf.h>
#include <typeinfo>

#include "readwrite.h"

using namespace netCDF;
// using namespace std;

void read_model(string fileModel, model *modP)
{
        cout << "reading model parameters from " << fileModel << "...";
        NcFile ParaInModel(fileModel.c_str(), NcFile::ReadOnly);

        NcVar* NpopP = ParaInModel.get_var("Npop");
        NpopP -> get(&modP->paras.Npop,1);

        modP->resize();
//         cout << "Npop: " << modP->paras.Npop << endl;
//         cout << "Npop type: " << typeid(modP->paras.Npop).name() << endl;
// get the membrane constants
        NcVar* tau_AP = ParaInModel.get_var("tau_A");
        NcVar* tau_NP = ParaInModel.get_var("tau_N");
//         NcVar* tau_GP = ParaInModel.get_var("tau_G");
        NcVar* tau_MP = ParaInModel.get_var("tau_M");

        tau_AP -> get(&modP->paras.tau_A,1);
        tau_NP -> get(&modP->paras.tau_N,1);
//         tau_GP -> get(&modP->paras.tau_G,1);
        tau_MP -> get(&modP->paras.tau_M,1);

// get the interaction parameters of populations
        NcVar* kappaP = ParaInModel.get_var("kappa");
//         NcVar* etaP = ParaInModel.get_var("eta");
//         NcVar* epsP = ParaInModel.get_var("eps");
//         NcVar* nP = ParaInModel.get_var("n");

        kappaP -> get(&modP->paras.kappa,1);
//         etaP -> get(&modP->paras.eta,1);
//         epsP -> get(&modP->paras.eps,1);
//         nP -> get(&modP->paras.n,1);

// get the drive parameters
        NcVar* driveP = ParaInModel.get_var("drive");
        NcVar* J_0P = ParaInModel.get_var("J_0");
        NcVar* K_0P = ParaInModel.get_var("K_0");
        NcVar* tau_0P = ParaInModel.get_var("tau_0");

        driveP -> get(&modP->paras.drive,1);
//         cout << "drive: " << modP->paras.drive << endl;
//         cout << "drive type: " << typeid(modP->paras.drive).name() << endl;

        if (modP->paras.drive > 0)    // if afferent, spiking drive, else constant input current
        {
                J_0P -> get(&modP->paras.J_0,1);
                K_0P -> get(&modP->paras.K_0,1);
                tau_0P -> get(&modP->paras.tau_0,1);
        }
        else
        {
                modP->paras.J_0 = 0;
                modP->paras.K_0 = 1;
                modP->paras.tau_0 = modP->paras.tau_A;
        }
//         cout << "read parameters from drive: rate=" << modP->paras.rateDrive << ", J_0=" << modP->paras.J_0 << ", K_0=" << modP->paras.K_0 << ", tau_0=" << modP->paras.tau_0 << endl;


// get simulation parameters
//         NcVar* rateWntP = ParaInModel.get_var("rateWnt");
//         NcVar* alpha_0P = ParaInModel.get_var("alpha_0");

//         rateWntP -> get(&modP->paras.rateWnt,1);
//         alpha_0P -> get(&modP->paras.alpha_0,1);
        cout << "done!" << endl;
}

void read_simulation(string fileSim, simulation *simP)
{
        cout << "reading simulation parameters from " << fileSim << "..." << endl;
        NcFile ParaInSim(fileSim.c_str(), NcFile::ReadOnly);

        NcDim* nSzP = ParaInSim.get_dim("nSz");
        NcDim* alpha_0SzP = ParaInSim.get_dim("alpha_0Sz");
        NcDim* tau_GSzP = ParaInSim.get_dim("tau_GSz");
        NcDim* rateWntSzP = ParaInSim.get_dim("rateWntSz");
        NcDim* epsSzP = ParaInSim.get_dim("epsSz");
        NcDim* etaSzP = ParaInSim.get_dim("etaSz");

        NcVar* nP = ParaInSim.get_var("n");
        NcVar* alpha_0P = ParaInSim.get_var("alpha_0");
        NcVar* tau_GP = ParaInSim.get_var("tau_G");
        NcVar* rateWntP = ParaInSim.get_var("rateWnt");
        NcVar* epsP = ParaInSim.get_var("eps");
        NcVar* etaP = ParaInSim.get_var("eta");

        NcVar* mode_calcP = ParaInSim.get_var("mode_calc");
        NcVar* mode_statsP = ParaInSim.get_var("mode_stats");

        simP->nSz = nSzP->size();
        simP->alpha_0Sz = alpha_0SzP->size();
        simP->tau_GSz = tau_GSzP->size();
        simP->rateWntSz = rateWntSzP->size();
        simP->epsSz = epsSzP->size();
        simP->etaSz = etaSzP->size();

//         cout << "sizes - n: " << simP->nSz << ", alpha_0: " << simP->alpha_0Sz << ", tau_G: " << simP->tau_GSz << ", rate: " << simP->rateWntSz << endl;

        mode_calcP -> get(&simP->mode_calc,1);
        mode_statsP -> get(&simP->mode_stats,1);

//         cout << "mode calc: " << simP->mode_calc << endl;
        simP->steps = max(max(max(simP->nSz,simP->alpha_0Sz),simP->tau_GSz),simP->rateWntSz);

//         cout << "steps: " << simP->steps << endl;

        simP -> n.resize(simP->nSz);
        simP -> alpha_0.resize(simP->alpha_0Sz);
        simP -> tau_G.resize(simP->tau_GSz);
        simP -> rateWnt.resize(simP->rateWntSz);
        simP -> eps.resize(simP->epsSz);
        simP -> eta.resize(simP->etaSz);

        nP -> get(&simP->n.front(),simP->nSz);
        if (simP->nSz == 1)
            simP->n.resize(simP->steps,simP->n[0]);

        alpha_0P -> get(&simP->alpha_0.front(),simP->alpha_0Sz);
        if (simP->alpha_0Sz == 1)
            simP -> alpha_0.resize(simP->steps,simP->alpha_0[0]);

        tau_GP -> get(&simP->tau_G.front(),simP->tau_GSz);
        if (simP->tau_GSz == 1)
            simP -> tau_G.resize(simP->steps,simP->tau_G[0]);

        rateWntP -> get(&simP->rateWnt.front(),simP->rateWntSz);
        epsP -> get(&simP->eps.front(),simP->epsSz);
        etaP -> get(&simP->eta.front(),simP->etaSz);
        cout << "done!" << endl;
}

void read_computation(string fileComp, computation *comP)
{

        cout << "reading computation parameters from " << fileComp << "..." ;
        NcFile ParaInComputation(fileComp.c_str(), NcFile::ReadOnly);

// get mode parameters
        NcVar* p_theoryP = ParaInComputation.get_var("p_theory");
        NcVar* p_theory_histP = ParaInComputation.get_var("p_theory_hist");
        NcVar* draw_from_theoryP = ParaInComputation.get_var("draw_from_theory");
        NcVar* draw_finite_timeP = ParaInComputation.get_var("draw_finite_time");
        NcVar* process_dataP = ParaInComputation.get_var("process_data");

        p_theoryP -> get(&comP->p_theory,1);
        p_theory_histP -> get(&comP->p_theory_hist,1);
        draw_from_theoryP -> get(&comP->draw_from_theory,1);
        draw_finite_timeP -> get(&comP->draw_finite_time,1);
        process_dataP -> get(&comP->process_data,1);

        NcVar* borderP = ParaInComputation.get_var("border");
        borderP -> get(&comP->border,1);

        NcVar* TP = ParaInComputation.get_var("T");
        TP -> get(&comP->T,1);
//         cout << "1" << endl;

        if (comP->draw_from_theory > 0)
        {
            NcDim* priorSzP = ParaInComputation.get_dim("priorSz");

            int priorSz = priorSzP->size();
            NcVar* priorP = ParaInComputation.get_var("prior");
            comP->prior.resize(priorSz);
            priorP -> get(&comP->prior.front(),priorSz);

            NcVar* NP = ParaInComputation.get_var("N");
            NP -> get(&comP->N,1);

            NcVar* n_binP = ParaInComputation.get_var("n_bin");
            n_binP -> get(&comP->n_bin,1);

            NcVar* seed_theoryP = ParaInComputation.get_var("seed_theory");
            comP -> seed_theory.resize(comP->draw_from_theory);
            seed_theoryP -> get(&comP->seed_theory.front(),comP->draw_from_theory);

            if (comP->draw_finite_time > 0)
            {
                NcVar* seed_timeP = ParaInComputation.get_var("seed_time");
                comP -> seed_time.resize(comP->draw_from_theory*comP->draw_finite_time);
                seed_timeP -> get(&comP->seed_time.front(),comP->draw_from_theory*comP->draw_finite_time);
            }
        }
        cout << "done!" << endl;
}



void read_measures(string fileMeasures, measures *mesP)
{
        cout << "reading measurement data from " << fileMeasures << "..." << endl;


        NcFile ParaInMeasures(fileMeasures.c_str(), NcFile::ReadOnly);

        NcDim* NSzP = ParaInMeasures.get_dim("NSz");
        mesP->N = NSzP->size();

        NcVar* TP = ParaInMeasures.get_var("T");
        TP -> get(&mesP->T,1);

//         int N_AP[NSz];
        mesP->N_AP.resize(mesP->N);
        NcVar* N_APP = ParaInMeasures.get_var("N_AP");
        N_APP -> get(&mesP->N_AP.front(),mesP->N);

        mesP->rates.resize(mesP->N);
        NcVar* ratesP = ParaInMeasures.get_var("rates");
        ratesP -> get(&mesP->rates.front(),mesP->N);

        cout << "done!" << endl;
}


/*
void write_results(string fileOut, simulation sim, double *gammaP, double *chiP)
{
	cout << "writing " << fileOut << "..." << endl;

	NcFile writeResults(fileOut.c_str(), NcFile::Replace);

	NcDim* two_dim = writeResults.add_dim("two", 2);

	NcVar *gamma = writeResults.add_var("gamma", ncDouble, two_dim);
	NcVar *chi = writeResults.add_var("chi", ncDouble, two_dim);

	gamma->put(gammaP, 2);
	chi->put(chiP, 2);
}*/

void write_theory(string fileOut, results *resP)
{
        cout << "(write theory) writing results to " << fileOut << "...";
        NcFile writeResults(fileOut.c_str(), NcFile::Replace);

        int p_Sz = resP->p_exact.size();
        NcDim* resolution_dim = writeResults.add_dim("resolution", p_Sz);

        NcVar *p_range = writeResults.add_var("p_range", ncDouble, resolution_dim);
        NcVar *p_exact = writeResults.add_var("p_exact", ncDouble, resolution_dim);
        NcVar *cdf_theory = writeResults.add_var("cdf_theory", ncDouble, resolution_dim);
//         NcVar *p_approx = writeResults.add_var("p_approx", ncDouble, resolution_dim);

        p_range->put(&resP->p_range.front(),p_Sz);
        p_exact->put(&resP->p_exact.front(),p_Sz);
        cdf_theory->put(&resP->cdf_theory.front(),p_Sz);
//         p_approx->put(&resP->p_approx.front(),p_Sz);

        cout << "done!" << endl;
}


// void write_theory(string fileOut, computation *comP, results *resP)
// {
//         cout << "writing results to " << fileOut << "...";
//
//         int p_Sz = resP->p_exact.size();
//         int bin_Sz = resP->p_hist.size();
//
//         NcFile writeResults(fileOut.c_str(), NcFile::Replace);
//
//         NcDim* one_dim = writeResults.add_dim("one", 1);
//
//         NcDim* k_dim = writeResults.add_dim("draw_from_theory", max(1,comP->draw_from_theory));
//         NcDim* j_dim = writeResults.add_dim("draw_finite_time", max(1,comP->draw_finite_time));
//
//         NcVar *d_nu = writeResults.add_var("d_nu", ncDouble, one_dim);
//         d_nu->put(&resP->d_nu,1);
//
//         if (comP->p_theory > 0)
//         {
//             NcDim* resolution_dim = writeResults.add_dim("resolution", p_Sz);
//
//             if (comP->p_theory >= 2)
//             {
//                 NcVar *alpha = writeResults.add_var("alpha", ncDouble, one_dim);
//                 NcVar *sigma_V = writeResults.add_var("sigma_V", ncDouble, one_dim);
//                 NcVar *I = writeResults.add_var("I", ncDouble, one_dim);
//                 NcVar *chi = writeResults.add_var("chi", ncDouble, one_dim);
//
//
//                 alpha->put(&resP->alpha_single,1);
//                 sigma_V->put(&resP->sigma_V_single,1);
//                 I->put(&resP->I_single,1);
//                 chi->put(&resP->chi_single,1);
//
//                 // write results from theory:
//                 NcVar *p_range = writeResults.add_var("p_range", ncDouble, resolution_dim);
//                 NcVar *p_exact = writeResults.add_var("p_exact", ncDouble, resolution_dim);
//                 NcVar *p_approx = writeResults.add_var("p_approx", ncDouble, resolution_dim);
//
//                 p_range->put(&resP->p_range.front(),p_Sz);
//                 p_exact->put(&resP->p_exact.front(),p_Sz);
//                 p_approx->put(&resP->p_approx.front(),p_Sz);
//
//                 if (comP->p_theory_hist > 0)
//                 {
//                     NcDim* bin_dim = writeResults.add_dim("bins", bin_Sz);
//
//                     NcVar *p_hist = writeResults.add_var("p_hist", ncDouble, bin_dim);
//                     p_hist->put(&resP->p_hist.front(),bin_Sz);
//                 }
//             }
//             // intersaving possible? to not crash netcdf files with too many datapoints
//
//     //         NcVar *p_k = writeResults.add_var("p_k", ncDouble, resolution_dim);
//     //         p_k->put(&resP->p_k.front(),p_Sz);
//
//             // write results from drawing rates
//
//
//             if (comP->draw_from_theory > 0)
//             {
//     //             NcVar *p_est_inf = writeResults.add_var("p_est_inf", ncDouble, resolution_dim);
//     //             p_est_inf->put(&resP->p_est_inf.front(),p_Sz);
//
//                 NcDim* N_dim = writeResults.add_dim("N", comP->N);
//
//                 NcVar *rate_inf = writeResults.add_var("rate_inf", ncDouble, k_dim, N_dim);
//
//                 double write_rate_inf[comP->draw_from_theory][comP->N];
//
//                 for (int k=0; k<comP->draw_from_theory; ++k)
//                     for (int n=0; n<comP->N; ++n)
//                     {
//     //                     cout << "rate_inf (" << k << "," << n << "): " << resP->rate_inf[k][n] << endl;
//                         write_rate_inf[k][n] = resP->rate_inf[k][n];
//                     }
//                 rate_inf->put(&write_rate_inf[0][0], comP->draw_from_theory, comP->N);
//             }
//
//             if (comP->draw_finite_time > 0)
//             {
//     //             NcVar *p_est_T = writeResults.add_var("p_est_T", ncDouble, resolution_dim);
//                 NcDim* N_dim = writeResults.add_dim("N", comP->N);
//
//                 NcVar *rate_T = writeResults.add_var("rate_T", ncDouble, k_dim, j_dim, N_dim);
//
//                 double write_rate_T[comP->draw_from_theory][comP->draw_finite_time][comP->N];
//                 for (int k=0; k<comP->draw_from_theory; ++k)
//                     for (int j=0; j<comP->draw_finite_time; ++j)
//                         for (int n=0; n<comP->N; ++n)
//                             write_rate_T[k][j][n] = resP->rate_T[k][j][n];
//
//                 rate_T->put(&write_rate_T[0][0][0], comP->draw_from_theory, comP->draw_finite_time, comP->N);
//
//
//                 NcVar *N_AP = writeResults.add_var("N_AP", ncDouble, k_dim, j_dim, N_dim);
//
//                 double write_N_AP[comP->draw_from_theory][comP->draw_finite_time][comP->N];
//                 for (int k=0; k<comP->draw_from_theory; ++k)
//                     for (int j=0; j<comP->draw_finite_time; ++j)
//                         for (int n=0; n<comP->N; ++n)
//                             write_N_AP[k][j][n] = resP->N_AP[k][j][n];
//
//                 N_AP->put(&write_N_AP[0][0][0], comP->draw_from_theory, comP->draw_finite_time, comP->N);
//
//
//
//                 NcVar *KS = writeResults.add_var("KS", ncDouble, k_dim, j_dim);
//
//                 double write_KS[comP->draw_from_theory][comP->draw_finite_time];
//                 for (int k=0; k<comP->draw_from_theory; ++k)
//                     for (int j=0; j<comP->draw_finite_time; ++j)
//                         write_KS[k][j] = resP->KS[k][j];
//
//                 KS->put(&write_KS[0][0], comP->draw_from_theory, comP->draw_finite_time);
//
//
//
// //                 NcVar *KL = writeResults.add_var("KL", ncDouble, k_dim, j_dim);
// //
// //                 double write_KL[comP->draw_from_theory][comP->draw_finite_time];
// //                 for (int k=0; k<comP->draw_from_theory; ++k)
// //                     for (int j=0; j<comP->draw_finite_time; ++j)
// //                         write_KL[k][j] = resP->KL_entropy[k][j];
//
// //                 KL->put(&write_KL[0][0], comP->draw_from_theory, comP->draw_finite_time);
//
//     //             p_est_T->put(&resP->p_est_T.front(),p_Sz);
//
//     //             double p_bayes_est_arr[comP->N][comP->AP_max];
//     //             for (int n=0; n<comP->N; ++n)
//     //             {
//     //                 for (int i=0; i<comP->AP_max; ++i)
//     //                     p_bayes_est_arr[n][i] = resP->p_bayes_est[n][i];
//     //             }
//
//     //             NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, N_dim, resolution_dim);
//
//                 NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, k_dim, j_dim, resolution_dim);
//
//                 double write_p_bayes_est[comP->draw_from_theory][comP->draw_finite_time][p_Sz];
//                 for (int k=0; k<comP->draw_from_theory; ++k)
//                     for (int j=0; j<comP->draw_finite_time; ++j)
//                         for (int i=0; i<p_Sz; ++i)
//                             write_p_bayes_est[k][j][i] = resP->p_bayes_est[k][j][i];
//
//                 p_bayes_est->put(&write_p_bayes_est[0][0][0], comP->draw_from_theory, comP->draw_finite_time, p_Sz);
//
//     //             NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, resolution_dim);
//     //             p_bayes_est->put(&resP->p_bayes_est.front(),p_Sz);
//     //             p_bayes_est->put(&p_bayes_est_arr[0][0], comP->N, p_Sz);
//             }
//         }
//         cout << "done!" << endl;
// }

void write_measures(string fileOut, computation *comP, measures *mesP, results *resP)
{
    //         cout << "writing rates to " << fileOut << "..." << endl;

        int p_Sz = resP->steps;
//         int bin_Sz = resP->p_hist.size();

        NcFile writeResults(fileOut.c_str(), NcFile::Replace);

        NcDim* one_dim = writeResults.add_dim("one", 1);
        NcDim* N_dim = writeResults.add_dim("N", mesP->N);
        NcDim* resolution_dim = writeResults.add_dim("resolution", p_Sz);
//         NcDim* bin_dim = writeResults.add_dim("bins", bin_Sz);

        NcVar *d_nu = writeResults.add_var("d_nu", ncDouble, one_dim);
        d_nu->put(&resP->d_nu,1);

//         // write results from theory:
//         NcVar *p_exact = writeResults.add_var("p_exact", ncDouble, resolution_dim);
//         NcVar *p_range = writeResults.add_var("p_range", ncDouble, resolution_dim);
//         NcVar *p_hist = writeResults.add_var("p_hist", ncDouble, bin_dim);
//
//         p_exact->put(&resP->p_exact.front(),p_Sz);
//         p_range->put(&resP->p_range.front(),p_Sz);
//         p_hist->put(&resP->p_hist.front(),bin_Sz);

        // intersaving possible? to not crash netcdf files with too many datapoints

        NcVar *rates = writeResults.add_var("rates", ncDouble, N_dim);
        rates->put(&mesP->rates.front(), mesP->N);

//         NcVar *p_k = writeResults.add_var("p_k", ncDouble, resolution_dim);
//         p_k->put(&resP->p_k.front(),p_Sz);
        NcVar *p_range = writeResults.add_var("p_range", ncDouble, resolution_dim);
        p_range->put(&resP->p_range.front(), p_Sz);

        NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, resolution_dim);
        p_bayes_est->put(&resP->p_bayes_est_measures.front(), p_Sz);

//             NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, resolution_dim);
//             p_bayes_est->put(&resP->p_bayes_est.front(),p_Sz);
//             p_bayes_est->put(&p_bayes_est_arr[0][0], comP->N, p_Sz);
}


void write_sharks(string fileOut, simulation sim, model mod, results *resP)//, double *gammaP, double *chiP, double *regionsP)
{
	cout << "writing simulation data to file '" << fileOut << "'..." << endl;

	NcFile writeResults(fileOut.c_str(), NcFile::Replace);

        int steps = resP->steps;
// 	cout << "size : " << resP->steps << endl;
// 	cout << "size: " << gammaP[0].size() << endl;
//         cout << "1" << endl;

// 	NcDim* Npop_dim = writeResults.add_dim("Npop", Npop);
	NcDim* steps_dim = writeResults.add_dim("steps", steps);

//         NcDim* rateWnt_dim = writeResults.add_dim("rateWntSz", sim.rateWntSz);
// 	NcDim* alpha_0_dim = writeResults.add_dim("alpha_0Sz", sim.alpha_0Sz);
// 	NcDim* tau_G_dim = writeResults.add_dim("tau_GSz", sim.tau_GSz);
// 	NcDim* n_dim = writeResults.add_dim("nSz", sim.nSz);
// 	NcDim* pop_dim = writeResults.add_dim("populations", mod.paras.population);

//         cout << "2" << endl;

	NcVar *gamma = writeResults.add_var("gamma", ncDouble,steps_dim, steps_dim);
	NcVar *chi = writeResults.add_var("chi", ncDouble, steps_dim, steps_dim);
        NcVar *regions = writeResults.add_var("regions", ncDouble, steps_dim, steps_dim);

        for (int rec=0; rec<resP->steps; ++rec)
        {
//                 cout << "plotting... " << endl;
//                 for (unsigned i=0; i<resP->steps; ++i)
//                     cout << "gamma: (" << rec << "," << i << "): " << resP->gamma[rec][i] << endl;
                gamma->put_rec(&resP->gamma[0][rec][0],rec);
                chi->put_rec(&resP->chi[0][rec][0],rec);
                regions->put_rec(&resP->regions[0][rec][0],rec);
        }

        if (mod.paras.Npop == 2)
        {
                NcVar *gamma_exc = writeResults.add_var("gamma_exc", ncDouble,steps_dim, steps_dim);
                NcVar *chi_exc = writeResults.add_var("chi_exc", ncDouble, steps_dim, steps_dim);
                NcVar *regions_exc = writeResults.add_var("regions_exc", ncDouble, steps_dim, steps_dim);

                for (int rec=0; rec<resP->steps; ++rec)
                {
        //                 cout << "plotting... " << endl;
        //                 for (unsigned i=0; i<resP->steps; ++i)
        //                     cout << "gamma: (" << rec << "," << i << "): " << resP->gamma[rec][i] << endl;
                        gamma_exc->put_rec(&resP->gamma[1][rec][0],rec);
                        chi_exc->put_rec(&resP->chi[1][rec][0],rec);
                        regions_exc->put_rec(&resP->regions[1][rec][0],rec);
                }
        }

// 	gamma->put(gammaP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);

// 	chi->put(chiP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);

//         regions->put(regionsP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);

//         cout << "3" << endl;
        if (sim.mode_stats == 4)
        {
                NcVar *gamma_approx = writeResults.add_var("gamma_approx", ncDouble, steps_dim, steps_dim);
                NcVar *chi_approx = writeResults.add_var("chi_approx", ncDouble, steps_dim, steps_dim);
                NcVar *entropy = writeResults.add_var("entropy", ncDouble, steps_dim, steps_dim);
                NcVar *KL = writeResults.add_var("KL", ncDouble, steps_dim, steps_dim);

                for (int rec=0; rec<resP->steps; ++rec)
                {
        //                 cout << "plotting... " << endl;
        //                 for (unsigned i=0; i<resP->steps; ++i)
        //                     cout << "gamma: (" << rec << "," << i << "): " << resP->gamma[rec][i] << endl;
                        gamma_approx->put_rec(&resP->gamma_approx[0][rec][0],rec);
                        chi_approx->put_rec(&resP->chi_approx[0][rec][0],rec);
//                         regions_approx->put_rec(&resP->regions_approx[rec][0],rec);
                        entropy->put_rec(&resP->entropy[0][rec][0],rec);
                        KL->put_rec(&resP->KL_entropy[0][rec][0],rec);
//                         regions_approx->put_rec(&resP->regions[rec][0],rec);
                }
//                     gamma_approx->put(gammaP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);

    //                 chi_approx->put(chiP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);

    //                 regions_approx->put(regionsP, sim.nSz, sim.alpha_0Sz, sim.tau_GSz, sim.rateWntSz, mod.paras.population);
        }


        NcVar *trans_DM = writeResults.add_var("DM_trans", ncDouble, steps_dim);
        NcVar *trans_np = writeResults.add_var("no_peak_trans", ncDouble, steps_dim);
        NcVar *trans_inc = writeResults.add_var("inc_trans", ncDouble, steps_dim);
        NcVar *trans_imp = writeResults.add_var("nu_implausible", ncDouble, steps_dim);

        int p_Sz = sim.alpha_0Sz;
//         cout << " steps: " << p_Sz << endl;

//                 cout << "DM1: " << resP->DM_trans[0] << " , DM90: " << resP->DM_trans[89] << endl;
        trans_DM->put(&resP->trans_DM[0][0], p_Sz);
        trans_np->put(&resP->trans_np[0][0], p_Sz);
        trans_inc->put(&resP->trans_inc[0][0], p_Sz);
        trans_imp->put(&resP->trans_imp[0][0], p_Sz);

        if (mod.paras.Npop == 2)
        {
                NcVar *trans_DM_exc = writeResults.add_var("DM_trans_exc", ncDouble, steps_dim);
                NcVar *trans_np_exc = writeResults.add_var("no_peak_trans_exc", ncDouble, steps_dim);
                NcVar *trans_inc_exc = writeResults.add_var("inc_trans_exc", ncDouble, steps_dim);
                NcVar *trans_imp_exc = writeResults.add_var("nu_implausible_exc", ncDouble, steps_dim);

                trans_DM_exc->put(&resP->trans_DM[1][0], p_Sz);
                trans_np_exc->put(&resP->trans_np[1][0], p_Sz);
                trans_inc_exc->put(&resP->trans_inc[1][0], p_Sz);
                trans_imp_exc->put(&resP->trans_imp[1][0], p_Sz);
        }

//         cout << "5" << endl;

//
// 	chi->put(resP->chi, resP->gamma.size(), resP->gamma[0].size(), 2);
}

void write_stats(string fileOut, model mod, simulation sim, results * resP)
{
        cout << "writing simulation data to file '" << fileOut << "'..." << endl;

//         static const int NC_ERR = 2;
//         NcError err(NcError::verbose_nonfatal);

	NcFile writeResults(fileOut.c_str(), NcFile::Replace);

        int dim1 = sim.max_ax[1];
        int dim2 = sim.max_ax[0];

        NcDim* Nc_dim1 = writeResults.add_dim("dim1", sim.max_ax[1]);
        NcDim* Nc_dim2 = writeResults.add_dim("dim2", sim.max_ax[0]);
//         NcDim* Nc_dim1 = writeResults.add_dim("dim1", sim.max_ax[1]);

        if (sim.mode_stats == 1)
        {
//                 write_stats_mode_1(&writeResults,resP,0,"",sim.max_ax[1],Nc_dim1,sim.max_ax[0],Nc_dim2);
//                 if sim.eps.size() > 1:
                NcVar *eps = writeResults.add_var("eps", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *eta = writeResults.add_var("eta", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *n = writeResults.add_var("n", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *tau_G = writeResults.add_var("tau_G", ncDouble, Nc_dim1, Nc_dim2);

//                 for (int i==0; i<dim2; ++i)
//                 {
                if (sim.eps.size() == 1)
                        sim.eps.resize(dim2,sim.eps[0]);
                if (sim.eta.size() == 1)
                        sim.eta.resize(dim2,sim.eta[0]);
                if (sim.n.size() == 1)
                        sim.n.resize(dim2,sim.n[0]);
                if (sim.tau_G.size() == 1)
                        sim.tau_G.resize(dim2,sim.tau_G[0]);


                NcVar *rateWnt = writeResults.add_var("rateWnt", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *q = writeResults.add_var("q", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *alpha_raw = writeResults.add_var("alpha_raw", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *alpha = writeResults.add_var("alpha", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *sigma_V = writeResults.add_var("sigma_V", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *gamma = writeResults.add_var("gamma", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *I_balance = writeResults.add_var("I_balance", ncDouble, Nc_dim1, Nc_dim2);
                NcVar *chi = writeResults.add_var("chi", ncDouble, Nc_dim1, Nc_dim2);

//                 cout << "now write" << endl;

                for (int rec=0; rec<dim1; ++rec)
                {
//                     cout << "eps size: " << sim.eps.size() << endl;
//                     cout << "eta size: " << sim.eta.size() << endl;
//                     cout << "n size: " << sim.n.size() << endl;
//                     cout << "tau_G size: " << sim.tau_G.size() << endl;

                    eps->put_rec(&sim.eps[0],rec);
                    eta->put_rec(&sim.eta[0],rec);
                    n->put_rec(&sim.n[0],rec);
                    tau_G->put_rec(&sim.tau_G[0],rec);

                    rateWnt->put_rec(&resP->rate[0][rec][0],rec);
                    q->put_rec(&resP->q[0][rec][0],rec);
                    alpha_raw->put_rec(&resP->alpha_raw[0][rec][0],rec);
                    alpha->put_rec(&resP->alpha[0][rec][0],rec);
                    sigma_V->put_rec(&resP->sigma_V[0][rec][0],rec);
                    gamma->put_rec(&resP->gamma[0][rec][0],rec);
                    I_balance->put_rec(&resP->I_balance[0][rec][0],rec);
                    chi->put_rec(&resP->chi[0][rec][0],rec);
                }
                NcVar *trans_imp = writeResults.add_var("nu_implausible", ncDouble, Nc_dim1);
                trans_imp->put(&resP->trans_imp[0][0], dim1);


                if (mod.paras.Npop == 2)
                {
                    cout << "#2" << endl;

                    NcVar *rateWnt = writeResults.add_var("rateWnt_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *q = writeResults.add_var("q_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *alpha_raw = writeResults.add_var("alpha_raw_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *alpha = writeResults.add_var("alpha_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *sigma_V = writeResults.add_var("sigma_V_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *gamma = writeResults.add_var("gamma_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *I_balance = writeResults.add_var("I_balance_exc", ncDouble, Nc_dim1, Nc_dim2);
                    NcVar *chi = writeResults.add_var("chi_exc", ncDouble, Nc_dim1, Nc_dim2);
                    cout << "now write" << endl;

                    for (int rec=0; rec<dim1; ++rec)
                    {
                        rateWnt->put_rec(&resP->rate[1][rec][0],rec);
                        q->put_rec(&resP->q[1][rec][0],rec);
                        alpha_raw->put_rec(&resP->alpha_raw[1][rec][0],rec);
                        alpha->put_rec(&resP->alpha[1][rec][0],rec);
                        sigma_V->put_rec(&resP->sigma_V[1][rec][0],rec);
                        gamma->put_rec(&resP->gamma[1][rec][0],rec);
                        I_balance->put_rec(&resP->I_balance[1][rec][0],rec);
                        chi->put_rec(&resP->chi[1][rec][0],rec);
                    }

                    NcVar *trans_imp = writeResults.add_var("nu_implausible_exc", ncDouble, Nc_dim1);
                    trans_imp->put(&resP->trans_imp[1][0], dim1);


//                     write_stats_mode_1(&writeResults,resP,1,"_exc",sim.max_ax[1],Nc_dim1,sim.max_ax[0],Nc_dim2);
                }
        }

        if (sim.mode_stats == 2)
        {
                NcDim* alpha_0_dim = writeResults.add_dim("alpha_0Sz", sim.alpha_0Sz);

                write_stats_mode_2(&writeResults,alpha_0_dim,resP,0,"");

                if (mod.paras.Npop == 2)
                        write_stats_mode_2(&writeResults,alpha_0_dim,resP,1,"_exc");
        }
}
//         if (sim.mode_stats == 3)
//         {
// //                 int p_Sz = resP->steps;
//                 NcDim* rateWnt_dim = writeResults.add_dim("rateWntSz", sim.rateWntSz);
//                 NcDim* alpha_0_dim = writeResults.add_dim("alpha_0Sz", sim.alpha_0Sz);
//
//
//                 NcVar *entropy = writeResults.add_var("entropy", ncDouble, alpha_0_dim,rateWnt_dim);
//                 NcVar *KL = writeResults.add_var("KL", ncDouble, alpha_0_dim,rateWnt_dim);
//
//                 for (unsigned rec=0; rec<sim.alpha_0Sz; ++rec)
//                 {
//                         entropy->put_rec(&resP->entropy[rec][0], rec);
//                         KL->put_rec(&resP->KL_entropy[rec][0], rec);
//                 }
//
//                 NcVar *rateWnt = writeResults.add_var("rateWnt", ncDouble, alpha_0_dim, rateWnt_dim);
//                 NcVar *q = writeResults.add_var("q", ncDouble, alpha_0_dim, rateWnt_dim);
// //                 NcVar *alpha = writeResults.add_var("alpha", ncDouble, alpha_0_dim, rateWnt_dim);
//                 NcVar *gamma = writeResults.add_var("gamma", ncDouble, alpha_0_dim, rateWnt_dim);
//                 NcVar *chi = writeResults.add_var("chi", ncDouble, alpha_0_dim, rateWnt_dim);
//
//                 NcVar *q_approx = writeResults.add_var("q_approx", ncDouble, alpha_0_dim, rateWnt_dim);
//                 NcVar *gamma_approx = writeResults.add_var("gamma_approx", ncDouble, alpha_0_dim, rateWnt_dim);
//                 NcVar *chi_approx = writeResults.add_var("chi_approx", ncDouble, alpha_0_dim, rateWnt_dim);
//                 for (unsigned rec=0; rec<sim.alpha_0Sz; ++rec)
//                 {
//                     rateWnt->put_rec(&resP->rateWnt[rec][0],rec);
//                     q->put_rec(&resP->q[rec][0],rec);
// //                     alpha->put_rec(&resP->alpha[rec][0],rec);
//                     gamma->put_rec(&resP->gamma[rec][0],rec);
//                     chi->put_rec(&resP->chi[rec][0],rec);
//
//                     q_approx->put_rec(&resP->q_approx[rec][0],rec);
//                     gamma_approx->put_rec(&resP->gamma_approx[rec][0],rec);
//                     chi_approx->put_rec(&resP->chi_approx[rec][0],rec);
//                 }
//         }
// }


// void write_stats_mode_1(NcFile *writeResults, results *resP, int p, string addStr, int dim1, NcDim* Nc_dim1, int dim2, NcDim* Nc_dim2)
// {
//         cout << "rate" << endl;
//
// //         if (p == 0)
// //         {
// //         NcDim* Nc_dim1 = writeResults->add_dim("alpha_0Sz", dim1);
// //         NcDim* Nc_dim2 = writeResults->add_dim("rateWntSz", dim2);
// //         }
//
// //         cout << "rateWnt"+addStr << endl;
// //         string str;
//
//         NcVar *rateWnt = writeResults->add_var("rateWnt", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *q = writeResults->add_var("q", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *alpha_raw = writeResults->add_var("alpha_raw", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *alpha = writeResults->add_var("alpha", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *sigma_V = writeResults->add_var("sigma_V", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *gamma = writeResults->add_var("gamma", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *I_balance = writeResults->add_var("I_balance", ncDouble, Nc_dim1, Nc_dim2);
//         NcVar *chi = writeResults->add_var("chi", ncDouble, Nc_dim1, Nc_dim2);
//         cout << "now write" << endl;
//
//         cout << "size rate: " << resP->rate.size() << "/" << resP->rate[0].size() << "/" << resP->rate[0][0].size() << endl;
//         cout << "size q: " << resP->q.size() << "/" << resP->q[0].size() << "/" << resP->q[0][0].size() << endl;
//         cout << "size alpha: " << resP->alpha.size() << "/" << resP->alpha[0].size() << "/" << resP->alpha[0][0].size() << endl;
//         for (int rec=0; rec<dim1; ++rec)
//         {
//             rateWnt->put_rec(&resP->rate[p][rec][0],rec);
//             q->put_rec(&resP->q[p][rec][0],rec);
//             alpha_raw->put_rec(&resP->alpha_raw[p][rec][0],rec);
//             alpha->put_rec(&resP->alpha[p][rec][0],rec);
//             sigma_V->put_rec(&resP->sigma_V[p][rec][0],rec);
//             gamma->put_rec(&resP->gamma[p][rec][0],rec);
//             I_balance->put_rec(&resP->I_balance[p][rec][0],rec);
//             chi->put_rec(&resP->chi[p][rec][0],rec);
//         }
//         cout << " done..." << endl;
//         NcVar *trans_imp = writeResults->add_var("nu_implausible", ncDouble, Nc_dim1);
//         trans_imp->put(&resP->trans_imp[p][0], dim1);
//         cout << " yo" << endl;
// }

void write_stats_mode_2(NcFile *writeResults, NcDim *Nc_dim, results *resP, int p, string addStr)
{
//         cout << "#2" << endl;
//     cout << strApp("inc_trans",addStr) << endl;
        NcVar *trans_DM = writeResults->add_var("DM_trans", ncDouble, Nc_dim);
        NcVar *trans_np = writeResults->add_var("no_peak_trans", ncDouble, Nc_dim);
        NcVar *trans_inc = writeResults->add_var("inc_trans", ncDouble, Nc_dim);
        NcVar *trans_imp = writeResults->add_var("nu_implausible", ncDouble, Nc_dim);

//         cout << resP->trans_DM.size() << endl;
//         cout << resP->trans_inc.size() << endl;
        trans_DM->put(&resP->trans_DM[p][0], resP->steps);
        trans_np->put(&resP->trans_np[p][0], resP->steps);
        trans_inc->put(&resP->trans_inc[p][0], resP->steps);
        trans_imp->put(&resP->trans_imp[p][0], resP->steps);
}

string strApp(string str1, string str2)
{
        string str="";
        str.append(str1);
        str.append(str2);
        return str;
}
