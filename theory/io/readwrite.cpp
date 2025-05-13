#include <iostream>
// #include <netcdfcpp.h>
#include <netcdf.h>
#include <typeinfo>

#include <plog/Log.h>

#include "readwrite.h"

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

bool get_from_ncid(int ncid, const char *varName, void *varP)
{
    // function to obtain information on variable values
    int varid;
    int status = nc_inq_varid(ncid, varName, &varid);
    if (status != NC_NOERR)
    {
        PLOG_INFO << "Variable "<< varName << " not found!";
        return false;
    }
    nc_get_var(ncid, varid, varP);
    return true;
};

bool get_from_ncid(int ncid, const char *varName, size_t *varP)
{
    // function to obtain information on variable dimension
    int varid;
    int status = nc_inq_dimid(ncid, varName, &varid);
    if (status != NC_NOERR)
    {
        PLOG_INFO << "Dimension "<< varName << " not found!";
        return false;
    }
    nc_inq_dimlen(ncid, varid, varP);
    return true;
};

template <typename T>
void get_dim_and_read(int ncid, const char *varName, size_t *varSz, vector<T> *varP)
{
    bool status = get_from_ncid(ncid, strApp(varName,"Sz").c_str(), varSz);
    if (status) {
        varP->resize(*varSz);
        status = get_from_ncid(ncid, varName, &varP->front());
    }
}

void write_to_ncid(int ncid, const char *varName, nc_type varType, int dimSz, const int *dimids, double *varData)
{
    int varid;
    nc_def_var(ncid, varName, varType, dimSz, dimids, &varid);
    nc_put_var(ncid, varid, varData);
};

void read_model(string fileModel, Model *modP)
{
    PLOG_INFO << "Reading 'model' parameters done from" << fileModel << "... ";

    int ncid;
    int status = true;

    nc_open(fileModel.c_str(), NC_NOWRITE, &ncid);

    status = status && get_from_ncid(ncid, "L", &modP->L);
    PLOG_DEBUG << "L: " << modP->L << endl;
    modP->layer.resize(modP->L);

    vector<unsigned> S;
    size_t nP;
    get_dim_and_read(ncid, "S", &nP, &S);

    unsigned nS = 0;
    for (unsigned p=0; p<nP; p++)
        nS += S[p];
    modP->nPop = nP;

    status = status && get_from_ncid(ncid, "I_alpha", &modP->I_alpha);
    status = status && get_from_ncid(ncid, "I_beta", &modP->I_beta);

    double eps[modP->L], eta[modP->L], J0_l[modP->L][modP->L];
    status = status && get_from_ncid(ncid, "eps", &eps[0]);
    status = status && get_from_ncid(ncid, "eta", &eta[0]);
    status = status && get_from_ncid(ncid, "J0_l", &J0_l[0]);

    double I_ext[nP];
    double rateWnt[nP], kappa[nP], alpha_0[nP], Psi_0[nP], tau_M[nP], tau_n[nP], J0[nP];
    status = status && get_from_ncid(ncid, "I_ext", &I_ext[0]);
    status = status && get_from_ncid(ncid, "rateWnt", &rateWnt[0]);
    status = status && get_from_ncid(ncid, "kappa", &kappa[0]);
    status = status && get_from_ncid(ncid, "alpha_0", &alpha_0[0]);
    status = status && get_from_ncid(ncid, "Psi_0", &Psi_0[0]);
    status = status && get_from_ncid(ncid, "tau_M", &tau_M[0]);
    status = status && get_from_ncid(ncid, "tau_n", &tau_n[0]);
    status = status && get_from_ncid(ncid, "J0", &J0[0]);

    double tau_I[nS], tau_norm[nS];
    status = status && get_from_ncid(ncid, "tau_I", &tau_I[0]);
    status = status && get_from_ncid(ncid, "tau_norm", &tau_norm[0]);

    unsigned p_idx = 0, s_idx = 0;
    for (unsigned l=0; l<modP->L; l++) {
        modP->layer[l].l = l;
        modP->layer[l].nPop = nP;
        modP->layer[l].population.resize(modP->layer[l].nPop);

        modP->layer[l].eps = eps[l];
        modP->layer[l].eta = eta[l];
        modP->layer[l].J0_l.resize(modP->L);
        for (unsigned ll=0; ll<modP->L; ll++)
            modP->layer[l].J0_l[ll] = J0_l[ll][l];    // check again, whether indices are in right order
        // modP->layer[l].setWeights();

        for (unsigned p=0; p<modP->layer[l].nPop; p++)
        {
            modP->layer[l].population[p].p = p;
            modP->layer[l].population[p].psp.resize(S[p_idx]);
            modP->layer[l].population[p].J.resize(modP->layer[l].nPop);

            modP->layer[l].population[p].nPSP = S[p_idx];

            modP->layer[l].population[p].I_ext = I_ext[p_idx];
            PLOG_DEBUG << "I_ext (" << l << "," << p << "): " << modP->layer[l].population[p].I_ext << " vs " << I_ext[p_idx] << endl;
            modP->layer[l].population[p].rateWnt = rateWnt[p_idx];
            modP->layer[l].population[p].kappa = kappa[p_idx];
            modP->layer[l].population[p].alpha_0 = alpha_0[p_idx];
            modP->layer[l].population[p].Psi_0 = Psi_0[p_idx];
            modP->layer[l].population[p].tau_M = tau_M[p_idx];
            modP->layer[l].population[p].tau_n = tau_n[p_idx];
            modP->layer[l].population[p].J0 = J0[p_idx];

            for (unsigned s=0; s<S[p_idx]; s++) {
                modP->layer[l].population[p].psp[s].s = s;
                modP->layer[l].population[p].psp[s].tau_I = tau_I[s_idx];
                modP->layer[l].population[p].psp[s].tau_norm = tau_norm[s_idx];

                s_idx++;
            }
            p_idx++;
        }
        // modP->layer[l].print_layer();
    }


// get the drive parameters
    // get_from_ncid(ncid, "drive", &modP->paras.drive);
    //
    // if (modP->paras.drive > 0)    // if afferent, spiking drive, else constant input current
    // {
    //     get_from_ncid(ncid, "J_0", &modP->paras.J_0);
    //     get_from_ncid(ncid, "K_0", &modP->paras.K_0);
    //     get_from_ncid(ncid, "tau_0", &modP->paras.tau_0);
    // } else {
    //     modP->paras.J_0 = 0;
    //     modP->paras.K_0 = 1;
    //     modP->paras.tau_0 = modP->paras.tau_A;
    // }

    nc_close(ncid);
    if (!status) throw;

    PLOG_INFO << "Reading 'model' done";
}

void read_simulation(string fileSim, Simulation *simP)
{
    PLOG_INFO << "Reading 'simulation' parameters done from" << fileSim << "... ";
    int ncid;
    int status = true;

    nc_open(fileSim.c_str(), NC_NOWRITE, &ncid);

    get_dim_and_read(ncid, "eps", &simP->epsSz, &simP->eps);
    get_dim_and_read(ncid, "eta", &simP->etaSz, &simP->eta);
    get_dim_and_read(ncid, "alpha_0", &simP->alpha_0Sz, &simP->alpha_0);
    get_dim_and_read(ncid, "Psi_0", &simP->Psi_0Sz, &simP->Psi_0);
    get_dim_and_read(ncid, "rateWnt", &simP->rateWntSz, &simP->rateWnt);
    get_dim_and_read(ncid, "tau_I", &simP->tau_ISz, &simP->tau_I);
    get_dim_and_read(ncid, "tau_n", &simP->tau_nSz, &simP->tau_n);
    get_dim_and_read(ncid, "I_alpha", &simP->I_alphaSz, &simP->I_alpha);
    get_dim_and_read(ncid, "I_beta", &simP->I_betaSz, &simP->I_beta);

    simP->sim_pointer.resize(2);
    get_dim_and_read(ncid, "sim_prim", &simP->sim_primSz, &simP->sim_pointer[0]);
    get_dim_and_read(ncid, "sim_sec", &simP->sim_secSz, &simP->sim_pointer[1]);
    
    status = status && get_from_ncid(ncid, "orderSz", &simP->orderSz);
    status = status && get_from_ncid(ncid, "charSz", &simP->charSz);

    PLOG_DEBUG << "orderSz: " << simP->orderSz << ", charSz: " << simP->charSz;

    simP->order.resize(simP->orderSz);

    // read variables from ncid
    status = status && get_from_ncid(ncid, "mode_calc", &simP->mode_calc);
    status = status && get_from_ncid(ncid, "mode_stats", &simP->mode_stats);
    status = status && get_from_ncid(ncid, "mode_selfcon", &simP->mode_selfcon);

    // reading string... quite bulky - isn't there a better solution?
    vector<vector<char> > order(simP->orderSz, vector<char> (simP->charSz));
    int varid;
    status = status && nc_inq_varid(ncid, "order", &varid);
    size_t start[2] = {0,0}, count[2] = {1,simP->charSz};

    PLOG_DEBUG << "Reading order variables... " << endl;
    for (unsigned rec=0;rec<simP->orderSz;rec++)//
    {
        start[0] = rec;
        nc_get_vara(ncid, varid, start, count, &order[rec][0]);
        simP->order[rec] = "";
        for (unsigned i=0;i<simP->charSz;i++) {
            simP->order[rec] += order[rec][i];
        }
        PLOG_DEBUG << simP->order[rec] << endl;
    }

    nc_close(ncid);
    PLOG_INFO << "Reading 'simulation' done";
}

void read_computation(string fileComp, Computation *comP)
{
    PLOG_INFO << "Reading 'computation' parameters done from" << fileComp << "... ";
    
    int ncid;
    nc_open(fileComp.c_str(), NC_NOWRITE, &ncid);

    // get mode parameters
    // get_from_ncid(ncid, "p_theory", &comP->p_theory);
    // get_from_ncid(ncid, "p_theory_hist", &comP->p_theory_hist);
    // get_from_ncid(ncid, "process_data", &comP->process_data);

    // get_from_ncid(ncid, "border", &comP->border);
    get_from_ncid(ncid, "draw_from_theory", &comP->draw_from_theory);
    get_from_ncid(ncid, "draw_finite_time", &comP->draw_finite_time);
    
    get_from_ncid(ncid, "N", &comP->N);
    get_from_ncid(ncid, "T", &comP->T);      // if T=0 => infty

    get_from_ncid(ncid, "seed", &comP->seed);
    // cout << "nTheory: " << comP->draw_from_theory << endl;
    // cout << "nTime: " << comP->draw_finite_time << endl;
    // cout << "seed: " << comP->seed << endl;
    // get_dim_and_read(ncid, "seed_theory", &comP->draw_from_theory, &comP->seed_theory);
    // get_dim_and_read(ncid, "seed_time", &comP->draw_finite_time, &comP->seed_time);
    // comP->seed_theory.resize(comP->draw_from_theory);
    // get_from_ncid(ncid, "seed_theory", &comP->seed_theory);
    // comP->seed_time.resize(comP->draw_finite_time);
    // get_from_ncid(ncid, "seed_time", &comP->seed_time);


    // if (comP->draw_from_theory > 0)
    // {
        // int priorSz;
        // get_from_ncid(ncid, "priorSz", &priorSz);
        // comP->prior.resize(priorSz);

        // get_from_ncid(ncid, "prior", &comP->prior);
        // get_from_ncid(ncid, "N", &comP->N);
        // get_from_ncid(ncid, "n_bin", &comP->n_bin);

        // get_from_ncid(ncid, "seed_theory", &comP->seed_theory);

        // if (comP->draw_finite_time > 0)
        // {
            // comP -> seed_time.resize(comP->draw_from_theory*comP->draw_finite_time);
            // get_from_ncid(ncid, "seed_time", &comP->seed_time);
        // }
    // }
    nc_close(ncid);
    PLOG_INFO << "Reading 'computation' done";
    // cout << "done!" << endl;
}



// void read_measures(string fileMeasures, Measures *mesP)
// {
//     // cout << "reading measurement data from " << fileMeasures << "... ";
//
//     int ncid;
//     nc_open(fileMeasures.c_str(), NC_NOWRITE, &ncid);
//
//     get_from_ncid(ncid, "NSz", &mesP->N);
//     get_from_ncid(ncid, "T", &mesP->T);
//
//     mesP->N_AP.resize(mesP->N);
//     get_from_ncid(ncid, "N_AP", &mesP->N_AP);
//     mesP->rates.resize(mesP->N);
//     get_from_ncid(ncid, "rates", &mesP->rates);
//     nc_close(ncid);
//     // cout << "done!" << endl;
// }



// void write_theory(string fileOut, results *resP)
// {
//     cout << "(write theory) writing results to " << fileOut << "... ";
//
//     int ncid, resolution_dim;
//     nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);
//
//     int p_Sz = resP->p_exact.size();
//     nc_def_dim(ncid, "resolution", p_Sz, &resolution_dim);
//
//     write_to_ncid(ncid,"p_range", NC_DOUBLE, 1, &resolution_dim, &resP->p_range.front());
//     write_to_ncid(ncid,"p_exact", NC_DOUBLE, 1, &resolution_dim, &resP->p_exact.front());
//     write_to_ncid(ncid,"p_approx", NC_DOUBLE, 1, &resolution_dim, &resP->p_approx.front());
//     write_to_ncid(ncid,"cdf_theory", NC_DOUBLE, 1, &resolution_dim, &resP->cdf_theory.front());
//
//     nc_close(ncid);
//     cout << "done!" << endl;
// }



void write_measures(string fileOut, Computation *comP, Model *modP, Measures *mesP)
{
    PLOG_INFO << "Writing measures (results) to" << fileOut << "... ";
    
    int ncid, N_dim, P_dim, rates_dim, rates_T_dim;
    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    unsigned nK = comP->draw_from_theory;
    unsigned nT = comP->draw_finite_time;

    // int p_Sz = resP->steps;//, bin_Sz = resP->p_hist.size();

    // nc_def_dim(ncid, "one", 1, &one_dim);
    nc_def_dim(ncid, "P", modP->layer[0].nPop, &P_dim);
    nc_def_dim(ncid, "N", comP->N, &N_dim);
    nc_def_dim(ncid, "rates_dim", nK, &rates_dim);
    nc_def_dim(ncid, "rates_T_dim", nT, &rates_T_dim);

    int dimids[4] = {P_dim, N_dim, rates_dim, rates_T_dim};

    int rates_id, rates_T_id;
    nc_def_var(ncid, "rates", NC_DOUBLE, 3, &dimids[0], &rates_id);
    nc_def_var(ncid, "rates_T", NC_DOUBLE, 4, &dimids[0], &rates_T_id);

    nc_enddef(ncid);

    size_t start[] = {0,0,0,0}, count[] = {1, 1, nK, nT};
        // nc_def_dim(ncid, "resolution", p_Sz, &resolution_dim);
    // // nc_def_dim(ncid, "bin_dim", bin_Sz, &bin_dim);
    //
    // write_to_ncid(ncid,"d_nu", NC_DOUBLE, 1, &one_dim, &resP->d_nu);

    for (unsigned p=0; p<modP->layer[0].nPop; p++) {

        for (unsigned n=0; n<comP->N; n++) {
            start[0] = p;
            start[1] = n;
            start[2] = 0;

            count[2] = nK;
            nc_put_vara(ncid, rates_id, start, count, &mesP->rates[p][n][0]);

            for (unsigned k=0; k<nK; k++) {
                start[2] = k;

                count[2] = 1;
                // count[3] = nT;
                nc_put_vara(ncid, rates_T_id, start, count, &mesP->rates_T[p][n][k][0]);
            }
        }
    }
    // write_to_ncid(ncid,"rates", NC_DOUBLE, 1, &N_dim, &mesP->rates.front());
    // rates->put(&mesP->rates.front(), mesP->N);
    // // write_to_ncid(ncid,"p_range", NC_DOUBLE, 1, &resolution_dim, &resP->p_range.front());
    // // p_range->put(&resP->p_range.front(), p_Sz);
    // write_to_ncid(ncid,"p_bayes_est", NC_DOUBLE, 1, &resolution_dim, &resP->p_bayes_est_measures.front());  // whats with measures vs ~measures?
    // // p_bayes_est->put(&resP->p_bayes_est_measures.front(), p_Sz);
    // // write_to_ncid(ncid,"p_k", NC_DOUBLE, 1, &resolution_dim, &resP->p_k.front());
    // // NcVar *p_k = writeResults.add_var("p_k", ncDouble, resolution_dim);
    // // p_k->put(&resP->p_k.front(),p_Sz);
    nc_close(ncid);
    PLOG_INFO << "done!" << endl;
}


void write_results(string fileOut, Simulation *simP, Model *modP, Model *modP_approx)
{
	// spdlog::info("writing shark data to file '{}'...",fileOut);

    // could change this to separate functions, which individually open file in write access and add data?!

	PLOG_INFO << "writing result data to file " << fileOut << "...";

    int ncid, steps_dim, steps_dim1, Npop_dim, trans_dim;//, info_dim;
    unsigned steps = simP->vars[0].steps;
    unsigned steps_1 = simP->vars[1].steps;
    PLOG_DEBUG << "steps: " << steps << ", " << steps_1;

    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    // cout << "Npop: " << modP->nPop << endl;
    nc_def_dim(ncid, "Npop", modP->nPop, &Npop_dim);
    nc_def_dim(ncid, "steps_dim", steps, &steps_dim);
    nc_def_dim(ncid, "steps_dim1", steps_1, &steps_dim1);
    nc_def_dim(ncid, "trans_dim", modP->simulation.nTrans, &trans_dim);

    int dimids[3] = {Npop_dim, steps_dim1, steps_dim};
    const int trans_dimids[3] = {Npop_dim, steps_dim1, trans_dim};

    // cout << "trans_dimids: " << modP->nPop << "," << steps_1 << "," << (int) modP->simulation.nTrans << endl;

    int para_dimID[7];
    write_prep_paras(ncid,&para_dimID[0],simP);

    int DM_id, np_id, inc_id, imp_id;
    int DM_approx_id, np_approx_id, inc_approx_id, imp_approx_id;
    int q_id, gamma_id, delta_id, rate_max_id, chi_id, I_balance_id, regions_id, regions_approx_id, implausible_id, implausible_approx_id;
    int info_id;
    nc_def_var(ncid, "q", NC_DOUBLE, 3, &dimids[0], &q_id);
    nc_def_var(ncid, "gamma", NC_DOUBLE, 3, &dimids[0], &gamma_id);
    nc_def_var(ncid, "delta", NC_DOUBLE, 3, &dimids[0], &delta_id);
    nc_def_var(ncid, "rate_max", NC_DOUBLE, 3, &dimids[0], &rate_max_id);
    nc_def_var(ncid, "chi", NC_DOUBLE, 3, &dimids[0], &chi_id);
    nc_def_var(ncid, "I_balance", NC_DOUBLE, 3, &dimids[0], &I_balance_id);
    if ((simP->mode_stats == 0) || (simP->mode_stats == 3))
    {
        nc_def_var(ncid, "regions", NC_DOUBLE, 3, &dimids[0], &regions_id);
        if (simP->mode_stats == 3)
            nc_def_var(ncid, "regions_approx", NC_DOUBLE, 3, &dimids[0], &regions_approx_id);
        
        nc_def_var(ncid, "implausible", NC_DOUBLE, 3, &dimids[0], &implausible_id);
        if (simP->mode_stats == 3)
            nc_def_var(ncid, "implausible_approx", NC_DOUBLE, 3, &dimids[0], &implausible_approx_id);
    }
    nc_def_var(ncid, "inc_trans", NC_INT, 2, &trans_dimids[1], &inc_id);
    nc_def_var(ncid, "imp_trans", NC_INT, 2, &trans_dimids[1], &imp_id);
    nc_def_var(ncid, "DM_trans", NC_INT, 3, &trans_dimids[0], &DM_id);
    nc_def_var(ncid, "np_trans", NC_INT, 3, &trans_dimids[0], &np_id);

    int alpha_raw_id, alpha_id, sigma_V_id;
    if (simP->mode_stats == 1)
    {
        nc_def_var(ncid, "regions", NC_DOUBLE, 3, &dimids[0], &regions_id);
        nc_def_var(ncid, "implausible", NC_DOUBLE, 3, &dimids[0], &implausible_id);
        nc_def_var(ncid, "alpha_raw", NC_DOUBLE, 3, &dimids[0], &alpha_raw_id);
        nc_def_var(ncid, "alpha", NC_DOUBLE, 3, &dimids[0], &alpha_id);
        nc_def_var(ncid, "sigma_V", NC_DOUBLE, 3, &dimids[0], &sigma_V_id);
    }

    int q_approx_id, gamma_approx_id, chi_approx_id, entropy_id, KL_entropy_id;
    if ((simP->mode_stats == 2) || (simP->mode_stats == 3))
    {
        nc_def_var(ncid, "q_approx", NC_DOUBLE, 3, &dimids[0], &q_approx_id);
        nc_def_var(ncid, "gamma_approx", NC_DOUBLE, 3, &dimids[0], &gamma_approx_id);
        nc_def_var(ncid, "chi_approx", NC_DOUBLE, 3, &dimids[0], &chi_approx_id);
        nc_def_var(ncid, "entropy", NC_DOUBLE, 3, &dimids[0], &entropy_id);
        nc_def_var(ncid, "KL_entropy", NC_DOUBLE, 3, &dimids[0], &KL_entropy_id);

        nc_def_var(ncid, "inc_trans_approx", NC_INT, 2, &trans_dimids[1], &inc_approx_id);
        nc_def_var(ncid, "imp_trans_approx", NC_INT, 2, &trans_dimids[1], &imp_approx_id);
        nc_def_var(ncid, "DM_trans_approx", NC_INT, 3, &trans_dimids[0], &DM_approx_id);
        nc_def_var(ncid, "np_trans_approx", NC_INT, 3, &trans_dimids[0], &np_approx_id);
    }

    if (simP->mode_stats == 4)
    {
        nc_def_var(ncid, "infoContent", NC_DOUBLE, 3, &dimids[0], &info_id);
    }

    // int p_range_id, p_exact_id, p_approx_id, cdf_theory_id;
    // if (simP->mode_stats == 3)
    // {
    //     unsigned p_Sz = resP->p_exact.size();
    //     nc_def_dim(ncid, "distribution_dim", p_Sz, &distr_dim);
    //
    //     dimids[3] = distr_dim;
    //
    //     start[3] = 0;
    //     count[3] = p_Sz;
    //
    //     nc_def_var(ncid, "p_range", NC_DOUBLE, 4, &dimids[0], &p_range_id);
    //     nc_def_var(ncid, "p_exact", NC_DOUBLE, 4, &dimids[0], &p_exact_id);
    //     nc_def_var(ncid, "p_approx", NC_DOUBLE, 4, &dimids[0], &p_approx_id);
    //     nc_def_var(ncid, "cdf_theory", NC_DOUBLE, 4, &dimids[0], &cdf_theory_id);
    // }
    // cout << "def done " << endl;
    nc_enddef(ncid);

    write_paras(ncid, para_dimID, simP);

    // prepare trans-vectors for writing:
    Population_Results *popResP, *popResP_approx;

    size_t start[] = {0,0,0}, count[] = {steps_1, steps_1, steps};
    // cout << "counts: " << count[0] << ", " << count[1] << ", " << count[2] << endl;
    count[0] = 1;
    for (unsigned s=0;s<steps_1;s++) {
        start[0] = s;
        start[1] = 0;

        count[1] = modP->results.trans_inc[s].size();
        if (count[1] > 0) {
            nc_put_vara(ncid, inc_id, start, count, &modP->results.trans_inc[s][0]);
            
            count[1] = modP->results.trans_imp[s].size();
            nc_put_vara(ncid, imp_id, start, count, &modP->results.trans_imp[s][0]);
        }
    }


    if ((simP->mode_stats == 2) || (simP->mode_stats == 3)) {
        for (unsigned s=0;s<steps_1;s++) {
            start[0] = s;
            start[1] = 0;
            count[1] = modP_approx->results.trans_inc[s].size();
            nc_put_vara(ncid, inc_approx_id, start, count, &modP_approx->results.trans_inc[s][0]);
            count[1] = modP_approx->results.trans_imp[s].size();
            nc_put_vara(ncid, imp_approx_id, start, count, &modP_approx->results.trans_imp[s][0]);
        }
    }

    count[1] = 1;
    unsigned p_idx = 0;
    for (unsigned l = 0; l < modP->L; l++) {
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {

            popResP = &modP->layer[l].population[p].results;
            popResP_approx = &modP_approx->layer[l].population[p].results;

            start[0] = p_idx;

            for (unsigned s=0;s<steps_1;s++) {
                start[1] = s;
                count[2] = popResP->trans_DM[s].size();
                nc_put_vara(ncid, DM_id, start, count, &popResP->trans_DM[s][0]);
                count[2] = popResP->trans_np[s].size();
                nc_put_vara(ncid, np_id, start, count, &popResP->trans_np[s][0]);
            }

            if ((simP->mode_stats == 2) || (simP->mode_stats == 3))
            {
                for (unsigned s=0;s<steps_1;s++) {
                    start[1] = s;
                    count[2] = popResP_approx->trans_DM[s].size();
                    nc_put_vara(ncid, DM_approx_id, start, count, &popResP_approx->trans_DM[s][0]);
                    count[2] = popResP_approx->trans_np[s].size();
                    nc_put_vara(ncid, np_approx_id, start, count, &popResP_approx->trans_np[s][0]);
                }
            }

            count[2] = steps;
            for (unsigned rec=0; rec<steps_1; rec++) {
                start[1] = rec;
                start[2] = 0;
                // count[2] = steps;
                nc_put_vara(ncid, q_id, start, count, &popResP->q[rec][0]);
                nc_put_vara(ncid, gamma_id, start, count, &popResP->gamma[rec][0]);
                nc_put_vara(ncid, delta_id, start, count, &popResP->delta[rec][0]);
                nc_put_vara(ncid, rate_max_id, start, count, &popResP->rate_max[rec][0]);
                nc_put_vara(ncid, chi_id, start, count, &popResP->chi[rec][0]);
                nc_put_vara(ncid, I_balance_id, start, count, &popResP->I_balance[rec][0]);

                if ((simP->mode_stats == 0) || (simP->mode_stats == 3))
                {
                    nc_put_vara(ncid, regions_id, start, count, &popResP->regions[rec][0]);
                    if (simP->mode_stats == 3)
                        nc_put_vara(ncid, regions_approx_id, start, count, &popResP_approx->regions[rec][0]);

                    nc_put_vara(ncid, implausible_id, start, count, &popResP->implausible[rec][0]);
                    if (simP->mode_stats == 3)
                        nc_put_vara(ncid, implausible_approx_id, start, count, &popResP_approx->implausible[rec][0]);
                }
                if (simP->mode_stats == 1)
                {
                    nc_put_vara(ncid, regions_id, start, count, &popResP->regions[rec][0]);
                    nc_put_vara(ncid, implausible_id, start, count, &popResP->implausible[rec][0]);
                    nc_put_vara(ncid, alpha_raw_id, start, count, &popResP->alpha_raw[rec][0]);
                    nc_put_vara(ncid, alpha_id, start, count, &popResP->alpha[rec][0]);
                    nc_put_vara(ncid, sigma_V_id, start, count, &popResP->sigma_V[rec][0]);
                }
                if ((simP->mode_stats == 2) || (simP->mode_stats == 3))
                {
                    //cout << "writing: l=" << l << ", p=" << p << endl;
                    nc_put_vara(ncid, q_approx_id, start, count, &popResP_approx->q[rec][0]);
                    nc_put_vara(ncid, gamma_approx_id, start, count, &popResP_approx->gamma[rec][0]);
                    nc_put_vara(ncid, chi_approx_id, start, count, &popResP_approx->chi[rec][0]);
                    nc_put_vara(ncid, entropy_id, start, count, &popResP->entropy[rec][0]);
                    nc_put_vara(ncid, KL_entropy_id, start, count, &popResP->KL_entropy[rec][0]);
                    //cout << "done!" << endl;
                }
                if (simP->mode_stats == 4)
                {

                    // count[2] = 1;
                    nc_put_vara(ncid, info_id, start, count, &popResP->infoContent[rec][0]);
                    // for (unsigned rec1=0; rec1<steps; rec1++)
                    // {
                    //     // cout << "start: " << start[0] << "," << start[1] << "," << start[2] << "," << start[3] << endl;
                    //     // cout << "count: " << count[0] << "," << count[1] << "," << count[2] << "," << count[3] << endl;
                    //     start[2] = rec1;
                    //     nc_put_vara(ncid, info_id, start, count, &resP->infoContent[p][rec][rec1][0]);
                    // }
                }
    //         // if (simP->mode_stats == 3)
    //         // {
    //         //     for (unsigned rec1=0; rec1<steps; rec1++)
    //         //     {
    //         //         nc_put_vara(ncid, p_range_id, start, count, &resP->p_range[p][rec][rec1][0]);
    //         //         nc_put_vara(ncid, p_exact_id, start, count, &resP->p_exact[p][rec][rec1][0]);
    //         //         nc_put_vara(ncid, p_approx_id, start, count, &resP->p_approx[p][rec][rec1][0]);
    //         //         nc_put_vara(ncid, cdf_theory_id, start, count, &resP->cdf_theory[p][rec][rec1][0]);
    //         //     }
            }
            p_idx++;
        }
    }

    nc_close(ncid);
    PLOG_INFO << "writing results done!";
}

void write_prep_paras(int ncid, int *dimID, Simulation *simP)
{
    int dim;
    for (unsigned i=0;i<simP->nVar;i++)
    {
        // cout << "preparing " << simP->vars[i].name << endl;
        nc_def_dim(ncid, strApp(simP->vars[i].name,"_dim").c_str(), simP->vars[i].steps, &dim);
        nc_def_var(ncid, simP->vars[i].name.c_str(), NC_DOUBLE, 1, &dim, dimID+i);
    }
}

void write_paras(int ncid, int *dimID, Simulation *simP)
{
    for (unsigned i=0;i<simP->nVar;i++)
    {
        // cout << "writing " << simP->vars[i].name << endl;
        nc_put_var(ncid, *(dimID+i), simP->vars[i].valP);
    }
}

string strApp(string str1, string str2)
{
        return str1.append(str2);
}
