#include <iostream>
#include <netcdfcpp.h>

#include "read_single.h"

void read_model(string fileModel, model *modP)
{
        cout << "reading mode parameters from " << fileModel << "...";
        NcFile ParaInModel(fileModel.c_str(), NcFile::ReadOnly);
        
// get the membrane constants
        NcVar* tau_AP = ParaInModel.get_var("tau_A");
        NcVar* tau_NP = ParaInModel.get_var("tau_N");
        NcVar* tau_GP = ParaInModel.get_var("tau_G");
        NcVar* tau_MP = ParaInModel.get_var("tau_M");
        
        tau_AP -> get(&modP->paras.tau_A,1);
        tau_NP -> get(&modP->paras.tau_N,1);
        tau_GP -> get(&modP->paras.tau_G,1);
        tau_MP -> get(&modP->paras.tau_M,1);
        
// get the interaction parameters of populations
        NcVar* kappaP = ParaInModel.get_var("kappa");
        NcVar* etaP = ParaInModel.get_var("eta");
        NcVar* epsP = ParaInModel.get_var("eps");
        NcVar* nP = ParaInModel.get_var("n");
        
        kappaP -> get(&modP->paras.kappa,1);
        etaP -> get(&modP->paras.eta,1);
        epsP -> get(&modP->paras.eps,1);
        nP -> get(&modP->paras.n,1);
        
        // get simulation parameters
        NcVar* rateWntP = ParaInModel.get_var("rateWnt");
        NcVar* alpha_0P = ParaInModel.get_var("alpha_0");
        
        rateWntP -> get(&modP->paras.rateWnt,1);
        alpha_0P -> get(&modP->paras.alpha_0,1);
        cout << "done!" << endl;
}

void read_simulation(string fileSim, simulation *simP)
{
        cout << "reading simulation parameters from " << fileModel << "...";
        NcFile ParaInSim(fileSim.c_str(), NcFile::ReadOnly);
        
        NcDim* nSzP = ParaInSim.get_dim("nSz");
        NcDim* alpha_0SzP = ParaInSim.get_dim("alpha_0Sz");
        NcDim* tau_GSzP = ParaInsim.get_dim("tau_GSz");
        NcDim* rateWntSzP = ParaInSim.get_dim("rateWntSz");
        
        NcVar* nP = ParaInSim.get_var("n");
        NcVar* alpha_0P = ParaInSim.get_var("alpha_0");
        NcVar* tau_GP = ParaInSim.get_var("tau_G");
        NcVar* rateWntP = ParaInSim.get_var("rateWnt");
        
        simP->nSz = nSzP.size();
        simP->alpha_0Sz = alpha_0SzP.size();
        simP->tau_GSz = tau_GSzP.size();
        simP->rateWntSz = rateWntSzP.size();
        
        nP -> get(&simP->n.front(),simP->nSz);
        alpha_0P -> get(&simP->alpha_0.front(),simP->alpha_0Sz);
        tau_GP -> get(&simP->tau_G.front(),simP->tau_GSz);
        rateWntP -> get(&simP->rateWnt.front(),simP->rateWntSz);
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
        
//         
        
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

void write_theory(string fileOut, computation *comP, results *resP)
{
        cout << "writing results to " << fileOut << "...";
            
        int p_Sz = resP->p_exact.size();
        int bin_Sz = resP->p_hist.size();
        
        NcFile writeResults(fileOut.c_str(), NcFile::Replace);
        
        NcDim* one_dim = writeResults.add_dim("one", 1);
        
        NcDim* k_dim = writeResults.add_dim("draw_from_theory", max(1,comP->draw_from_theory));
        NcDim* j_dim = writeResults.add_dim("draw_finite_time", max(1,comP->draw_finite_time));
        
        NcVar *d_nu = writeResults.add_var("d_nu", ncDouble, one_dim);
        d_nu->put(&resP->d_nu,1);
        
        if (comP->p_theory > 0)
        {
            NcDim* resolution_dim = writeResults.add_dim("resolution", p_Sz);
            
            if (comP->p_theory >= 2)
            {
                // write results from theory:
                NcVar *p_exact = writeResults.add_var("p_exact", ncDouble, resolution_dim);
                NcVar *p_range = writeResults.add_var("p_range", ncDouble, resolution_dim);

                p_exact->put(&resP->p_exact.front(),p_Sz);
                p_range->put(&resP->p_range.front(),p_Sz);

                if (comP->p_theory_hist > 0)
                {
                    NcDim* bin_dim = writeResults.add_dim("bins", bin_Sz);
                
                    NcVar *p_hist = writeResults.add_var("p_hist", ncDouble, bin_dim);
                    p_hist->put(&resP->p_hist.front(),bin_Sz);
                }
            }
            // intersaving possible? to not crash netcdf files with too many datapoints
            
    //         NcVar *p_k = writeResults.add_var("p_k", ncDouble, resolution_dim);
    //         p_k->put(&resP->p_k.front(),p_Sz);
            
            // write results from drawing rates
            
            
            if (comP->draw_from_theory > 0)
            {
    //             NcVar *p_est_inf = writeResults.add_var("p_est_inf", ncDouble, resolution_dim);
    //             p_est_inf->put(&resP->p_est_inf.front(),p_Sz);
                
                NcDim* N_dim = writeResults.add_dim("N", comP->N);
                
                NcVar *rate_inf = writeResults.add_var("rate_inf", ncDouble, k_dim, N_dim);
                
                double write_rate_inf[comP->draw_from_theory][comP->N];
                
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int n=0; n<comP->N; ++n)
                    {
    //                     cout << "rate_inf (" << k << "," << n << "): " << resP->rate_inf[k][n] << endl;
                        write_rate_inf[k][n] = resP->rate_inf[k][n];
                    }
                rate_inf->put(&write_rate_inf[0][0], comP->draw_from_theory, comP->N);
            }
            
            if (comP->draw_finite_time > 0)
            {
    //             NcVar *p_est_T = writeResults.add_var("p_est_T", ncDouble, resolution_dim);
                NcDim* N_dim = writeResults.add_dim("N", comP->N);
                
                NcVar *rate_T = writeResults.add_var("rate_T", ncDouble, k_dim, j_dim, N_dim);
                
                double write_rate_T[comP->draw_from_theory][comP->draw_finite_time][comP->N];
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int j=0; j<comP->draw_finite_time; ++j)
                        for (int n=0; n<comP->N; ++n)
                            write_rate_T[k][j][n] = resP->rate_T[k][j][n];
                
                rate_T->put(&write_rate_T[0][0][0], comP->draw_from_theory, comP->draw_finite_time, comP->N);
                
                
                NcVar *N_AP = writeResults.add_var("N_AP", ncDouble, k_dim, j_dim, N_dim);
                
                double write_N_AP[comP->draw_from_theory][comP->draw_finite_time][comP->N];
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int j=0; j<comP->draw_finite_time; ++j)
                        for (int n=0; n<comP->N; ++n)
                            write_N_AP[k][j][n] = resP->N_AP[k][j][n];
                
                N_AP->put(&write_N_AP[0][0][0], comP->draw_from_theory, comP->draw_finite_time, comP->N);
                
                
                
                NcVar *KS = writeResults.add_var("KS", ncDouble, k_dim, j_dim);
                
                double write_KS[comP->draw_from_theory][comP->draw_finite_time];
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int j=0; j<comP->draw_finite_time; ++j)
                        write_KS[k][j] = resP->KS[k][j];
                
                KS->put(&write_KS[0][0], comP->draw_from_theory, comP->draw_finite_time);
                
                
                
                NcVar *KL = writeResults.add_var("KL", ncDouble, k_dim, j_dim);
                
                double write_KL[comP->draw_from_theory][comP->draw_finite_time];
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int j=0; j<comP->draw_finite_time; ++j)
                        write_KL[k][j] = resP->KL_entropy[k][j];
                
                KL->put(&write_KL[0][0], comP->draw_from_theory, comP->draw_finite_time);
                
    //             p_est_T->put(&resP->p_est_T.front(),p_Sz);
                
    //             double p_bayes_est_arr[comP->N][comP->AP_max];
    //             for (int n=0; n<comP->N; ++n)
    //             {
    //                 for (int i=0; i<comP->AP_max; ++i)
    //                     p_bayes_est_arr[n][i] = resP->p_bayes_est[n][i];
    //             }
                        
    //             NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, N_dim, resolution_dim);
                
                NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, k_dim, j_dim, resolution_dim);
                
                double write_p_bayes_est[comP->draw_from_theory][comP->draw_finite_time][p_Sz];
                for (int k=0; k<comP->draw_from_theory; ++k)
                    for (int j=0; j<comP->draw_finite_time; ++j)
                        for (int i=0; i<p_Sz; ++i)
                            write_p_bayes_est[k][j][i] = resP->p_bayes_est[k][j][i];
                
                p_bayes_est->put(&write_p_bayes_est[0][0][0], comP->draw_from_theory, comP->draw_finite_time, p_Sz);
                
    //             NcVar *p_bayes_est = writeResults.add_var("p_bayes_est", ncDouble, resolution_dim);
    //             p_bayes_est->put(&resP->p_bayes_est.front(),p_Sz);
    //             p_bayes_est->put(&p_bayes_est_arr[0][0], comP->N, p_Sz);
            }
        }
        cout << "done!" << endl;
}

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