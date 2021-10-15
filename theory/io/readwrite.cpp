#include <iostream>
// #include <netcdfcpp.h>
#include <netcdf.h>
#include <typeinfo>

#include "readwrite.h"

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void get_from_ncid(int ncid, const char *varName, void *varP)
{
    int varid;
    nc_inq_varid(ncid, varName, &varid);
    nc_get_var(ncid, varid, varP);
};

void get_from_ncid(int ncid, const char *varName, size_t *varP)
{
    int varid;
    nc_inq_dimid(ncid, varName, &varid);
    nc_inq_dimlen(ncid, varid, varP);
};

void write_to_ncid(int ncid, const char *varName, nc_type varType, int dimSz, const int *dimids, double *varData)
{
    int varid;
    nc_def_var(ncid, varName, varType, dimSz, dimids, &varid);
    nc_put_var(ncid, varid, varData);
};

void read_model(string fileModel, model *modP)
{
    // spdlog::info("reading model parameters from {} ... ",fileModel);
    cout << "reading model parameters from " << fileModel << "... ";
    int ncid;

    nc_open(fileModel.c_str(), NC_NOWRITE, &ncid);

    get_from_ncid(ncid, "Npop", &modP->paras.Npop);
    modP->resize();

    // get the membrane constants
    get_from_ncid(ncid, "tau_A", &modP->paras.tau_A);
    get_from_ncid(ncid, "tau_N", &modP->paras.tau_N);
    get_from_ncid(ncid, "tau_M", &modP->paras.tau_M);
    get_from_ncid(ncid, "J", &modP->paras.J);
    cout << "J: " << modP->paras.J << endl;
    // get_from_ncid(ncid, "tau_G", &modP->paras.tau_G);

// get the interaction parameters of populations
    get_from_ncid(ncid, "kappa", &modP->paras.kappa);

// get the drive parameters
    get_from_ncid(ncid, "drive", &modP->paras.drive);

    if (modP->paras.drive > 0)    // if afferent, spiking drive, else constant input current
    {
        get_from_ncid(ncid, "J_0", &modP->paras.J_0);
        get_from_ncid(ncid, "K_0", &modP->paras.K_0);
        get_from_ncid(ncid, "tau_0", &modP->paras.tau_0);
    } else {
        modP->paras.J_0 = 0;
        modP->paras.K_0 = 1;
        modP->paras.tau_0 = modP->paras.tau_A;
    }

    nc_close(ncid);
    cout << "done!" << endl;
}

void read_simulation(string fileSim, simulation *simP)
{
    cout << "reading simulation parameters from " << fileSim << "... ";

    int ncid;

    nc_open(fileSim.c_str(), NC_NOWRITE, &ncid);

    // find size of dimensions
    get_from_ncid(ncid, "nSz", &simP->nSz);
    get_from_ncid(ncid, "alpha_0Sz", &simP->alpha_0Sz);
    get_from_ncid(ncid, "tau_GSz", &simP->tau_GSz);
    get_from_ncid(ncid, "rateWntSz", &simP->rateWntSz);
    get_from_ncid(ncid, "epsSz", &simP->epsSz);
    get_from_ncid(ncid, "etaSz", &simP->etaSz);
    get_from_ncid(ncid, "orderSz", &simP->orderSz);
    get_from_ncid(ncid, "charSz", &simP->charSz);

    simP->steps = max(max(max(simP->nSz,simP->alpha_0Sz),simP->tau_GSz),simP->rateWntSz);
    // cout << "sizes - n: " << simP->nSz << ", alpha_0: " << simP->alpha_0Sz << ", tau_G: " << simP->tau_GSz << ", rate: " << simP->rateWntSz << endl;
    // cout << "steps: " << simP->steps << endl;

    // adjust sizes of simP accordingly
    simP->n.resize(simP->nSz);
    simP->alpha_0.resize(simP->alpha_0Sz);
    simP->tau_G.resize(simP->tau_GSz);
    simP->rateWnt.resize(simP->rateWntSz);
    simP->eps.resize(simP->epsSz);
    simP->eta.resize(simP->etaSz);
    simP->order.resize(simP->orderSz);

    // read variables from ncid
    get_from_ncid(ncid, "mode_calc", &simP->mode_calc);
    get_from_ncid(ncid, "mode_stats", &simP->mode_stats);

    // reading string... quite bulky - isn't there a better solution?
    vector<vector<char>> order(simP->orderSz, vector<char> (simP->charSz));
    int varid;
    nc_inq_varid(ncid, "order", &varid);
    size_t start[2] = {0,0}, count[2] = {1,simP->charSz};

    for (unsigned rec=0;rec<simP->orderSz;rec++)//
    {
        start[0] = rec;
        nc_get_vara(ncid, varid, start, count, &order[rec][0]);
        simP->order[rec] = "";
        for (unsigned i=0;i<simP->charSz;i++)
            simP->order[rec] += order[rec][i];
    }

    get_from_ncid(ncid, "n", &simP->n[0]);
    get_from_ncid(ncid, "alpha_0", &simP->alpha_0[0]);
    get_from_ncid(ncid, "tau_G", &simP->tau_G[0]);
    get_from_ncid(ncid, "rateWnt", &simP->rateWnt[0]);
    get_from_ncid(ncid, "eps", &simP->eps[0]);
    get_from_ncid(ncid, "eta", &simP->eta[0]);

    nc_close(ncid);
    cout << "done!" << endl;
}

void read_computation(string fileComp, computation *comP)
{
    cout << "reading computation parameters from " << fileComp << "... " ;

    int ncid;
    nc_open(fileComp.c_str(), NC_NOWRITE, &ncid);

// get mode parameters
    get_from_ncid(ncid, "p_theory", &comP->p_theory);
    get_from_ncid(ncid, "p_theory_hist", &comP->p_theory_hist);
    get_from_ncid(ncid, "draw_from_theory", &comP->draw_from_theory);
    get_from_ncid(ncid, "draw_finite_time", &comP->draw_finite_time);
    get_from_ncid(ncid, "process_data", &comP->process_data);

    get_from_ncid(ncid, "border", &comP->border);
    get_from_ncid(ncid, "T", &comP->T);

    if (comP->draw_from_theory > 0)
    {
        int priorSz;
        get_from_ncid(ncid, "priorSz", &priorSz);
        comP->prior.resize(priorSz);

        get_from_ncid(ncid, "prior", &comP->prior);
        get_from_ncid(ncid, "N", &comP->N);
        get_from_ncid(ncid, "n_bin", &comP->n_bin);

        comP -> seed_theory.resize(comP->draw_from_theory);
        get_from_ncid(ncid, "seed_theory", &comP->seed_theory);

        if (comP->draw_finite_time > 0)
        {
            comP -> seed_time.resize(comP->draw_from_theory*comP->draw_finite_time);
            get_from_ncid(ncid, "seed_time", &comP->seed_time);
        }
    }
    nc_close(ncid);
    cout << "done!" << endl;
}



void read_measures(string fileMeasures, measures *mesP)
{
    cout << "reading measurement data from " << fileMeasures << "... ";

    int ncid;
    nc_open(fileMeasures.c_str(), NC_NOWRITE, &ncid);

    get_from_ncid(ncid, "NSz", &mesP->N);
    get_from_ncid(ncid, "T", &mesP->T);

    mesP->N_AP.resize(mesP->N);
    get_from_ncid(ncid, "N_AP", &mesP->N_AP);
    mesP->rates.resize(mesP->N);
    get_from_ncid(ncid, "rates", &mesP->rates);
    nc_close(ncid);
    cout << "done!" << endl;
}



void write_theory(string fileOut, results *resP)
{
    cout << "(write theory) writing results to " << fileOut << "... ";

    int ncid, resolution_dim;
    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    int p_Sz = resP->p_exact.size();
    nc_def_dim(ncid, "resolution", p_Sz, &resolution_dim);

    write_to_ncid(ncid,"p_range", NC_DOUBLE, 1, &resolution_dim, &resP->p_range.front());
    write_to_ncid(ncid,"p_exact", NC_DOUBLE, 1, &resolution_dim, &resP->p_exact.front());
    write_to_ncid(ncid,"p_approx", NC_DOUBLE, 1, &resolution_dim, &resP->p_approx.front());
    write_to_ncid(ncid,"cdf_theory", NC_DOUBLE, 1, &resolution_dim, &resP->cdf_theory.front());

    nc_close(ncid);
    cout << "done!" << endl;
}



void write_measures(string fileOut, computation *comP, measures *mesP, results *resP)
{
    cout << "writing measures (results) to " << fileOut << "..." << endl;

    int ncid, one_dim, N_dim, resolution_dim;
    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    int p_Sz = resP->steps;//, bin_Sz = resP->p_hist.size();

    nc_def_dim(ncid, "one", 1, &one_dim);
    nc_def_dim(ncid, "N", mesP->N, &N_dim);
    nc_def_dim(ncid, "resolution", p_Sz, &resolution_dim);
    // nc_def_dim(ncid, "bin_dim", bin_Sz, &bin_dim);

    write_to_ncid(ncid,"d_nu", NC_DOUBLE, 1, &one_dim, &resP->d_nu);
    write_to_ncid(ncid,"rates", NC_DOUBLE, 1, &N_dim, &mesP->rates.front());
    // rates->put(&mesP->rates.front(), mesP->N);
    write_to_ncid(ncid,"p_range", NC_DOUBLE, 1, &resolution_dim, &resP->p_range.front());
    // p_range->put(&resP->p_range.front(), p_Sz);
    write_to_ncid(ncid,"p_bayes_est", NC_DOUBLE, 1, &resolution_dim, &resP->p_bayes_est_measures.front());  // whats with measures vs ~measures?
    // p_bayes_est->put(&resP->p_bayes_est_measures.front(), p_Sz);
    // write_to_ncid(ncid,"p_k", NC_DOUBLE, 1, &resolution_dim, &resP->p_k.front());
    // NcVar *p_k = writeResults.add_var("p_k", ncDouble, resolution_dim);
    // p_k->put(&resP->p_k.front(),p_Sz);
    nc_close(ncid);
    cout << "done!" << endl;
}


void write_results(string fileOut, simulation *simP, model *modP, results *resP)
{
	// spdlog::info("writing shark data to file '{}'...",fileOut);
	cout << "writing result data to file " << fileOut << "...";

    int ncid, dimids[4], steps_dim, steps_dim1, Npop_dim, distr_dim;
    int steps_1 = simP->vars[1].steps;
    int steps = simP->vars[0].steps;

    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    nc_def_dim(ncid, "Npop", modP->paras.Npop, &Npop_dim);
    nc_def_dim(ncid, "steps_dim1", steps_1, &steps_dim1);
    nc_def_dim(ncid, "steps_dim", steps, &steps_dim);

    dimids[0] = Npop_dim;
    dimids[1] = steps_dim1;
    dimids[2] = steps_dim;

    size_t start[4], count[4];
    start[2] = 0;
    count[0] = 1;
    count[2] = steps;

    int para_dimID[6];
    write_prep_paras(ncid,&para_dimID[0],simP);

    int DM_id, np_id, inc_id, imp_id;
    int q_id, gamma_id, chi_id, regions_id;
    nc_def_var(ncid, "q", NC_DOUBLE, 3, &dimids[0], &q_id);
    nc_def_var(ncid, "gamma", NC_DOUBLE, 3, &dimids[0], &gamma_id);
    nc_def_var(ncid, "chi", NC_DOUBLE, 3, &dimids[0], &chi_id);
    if (simP->mode_stats == 0)
        nc_def_var(ncid, "regions", NC_DOUBLE, 3, &dimids[0], &regions_id);
    nc_def_var(ncid, "DM_trans", NC_DOUBLE, 2, &dimids[0], &DM_id);
    nc_def_var(ncid, "np_trans", NC_DOUBLE, 2, &dimids[0], &np_id);
    nc_def_var(ncid, "inc_trans", NC_DOUBLE, 2, &dimids[0], &inc_id);
    nc_def_var(ncid, "imp_trans", NC_DOUBLE, 2, &dimids[0], &imp_id);

    int alpha_raw_id, alpha_id, sigma_V_id, I_balance_id;
    if (simP->mode_stats == 1)
    {
        nc_def_var(ncid, "alpha_raw", NC_DOUBLE, 3, &dimids[0], &alpha_raw_id);
        nc_def_var(ncid, "alpha", NC_DOUBLE, 3, &dimids[0], &alpha_id);
        nc_def_var(ncid, "sigma_V", NC_DOUBLE, 3, &dimids[0], &sigma_V_id);
        nc_def_var(ncid, "I_balance", NC_DOUBLE, 3, &dimids[0], &I_balance_id);
    }

    int q_approx_id, gamma_approx_id, chi_approx_id, entropy_id, KL_entropy_id;
    if (simP->mode_stats == 2)
    {
        nc_def_var(ncid, "q_approx", NC_DOUBLE, 3, &dimids[0], &q_approx_id);
        nc_def_var(ncid, "gamma_approx", NC_DOUBLE, 3, &dimids[0], &gamma_approx_id);
        nc_def_var(ncid, "chi_approx", NC_DOUBLE, 3, &dimids[0], &chi_approx_id);
        nc_def_var(ncid, "entropy", NC_DOUBLE, 3, &dimids[0], &entropy_id);
        nc_def_var(ncid, "KL_entropy", NC_DOUBLE, 3, &dimids[0], &KL_entropy_id);
    }

    if (simP->mode_stats == 3)
    {
        int p_Sz = resP->p_exact.size();
        nc_def_dim(ncid, "distribution_dim", p_Sz, &distr_dim);

        dimids[3] = distr_dim;

        start[3] = 0;
        count[3] = p_Sz;

        nc_def_var(ncid, "p_range", NC_DOUBLE, 4, &dimids[0], &p_range_id);
        nc_def_var(ncid, "p_exact", NC_DOUBLE, 4, &dimids[0], &p_exact_id);
        nc_def_var(ncid, "p_approx", NC_DOUBLE, 4, &dimids[0], &p_approx_id);
        nc_def_var(ncid, "cdf_theory", NC_DOUBLE, 4, &dimids[0], &cdf_theory_id);
    }

    nc_enddef(ncid);

    write_paras(ncid, para_dimID, simP);


    for (int p=0; p<modP->paras.Npop; p++)
    {
        start[0] = p;

        start[1] = 0;
        count[1] = steps_1;
        nc_put_vara(ncid, DM_id, start, count, &resP->trans_DM[p][0]);
        nc_put_vara(ncid, np_id, start, count, &resP->trans_np[p][0]);
        nc_put_vara(ncid, inc_id, start, count, &resP->trans_inc[p][0]);
        nc_put_vara(ncid, imp_id, start, count, &resP->trans_imp[p][0]);

        count[1] = 1;

        for (int rec=0; rec<steps_1; rec++)
        {
            start[1] = rec;
            nc_put_vara(ncid, q_id, start, count, &resP->q[p][rec][0]);
            nc_put_vara(ncid, gamma_id, start, count, &resP->gamma[p][rec][0]);
            nc_put_vara(ncid, chi_id, start, count, &resP->chi[p][rec][0]);

            if (simP->mode_stats == 0)
                nc_put_vara(ncid, regions_id, start, count, &resP->regions[p][rec][0]);
            if (simP->mode_stats == 1)
            {
                nc_put_vara(ncid, alpha_raw_id, start, count, &resP->alpha_raw[p][rec][0]);
                nc_put_vara(ncid, alpha_id, start, count, &resP->alpha[p][rec][0]);
                nc_put_vara(ncid, sigma_V_id, start, count, &resP->sigma_V[p][rec][0]);
                nc_put_vara(ncid, I_balance_id, start, count, &resP->I_balance[p][rec][0]);
            }
            if (simP->mode_stats == 2)
            {
                nc_put_vara(ncid, q_approx_id, start, count, &resP->q_approx[p][rec][0]);
                nc_put_vara(ncid, gamma_approx_id, start, count, &resP->gamma_approx[p][rec][0]);
                nc_put_vara(ncid, chi_approx_id, start, count, &resP->chi_approx[p][rec][0]);
                nc_put_vara(ncid, entropy_id, start, count, &resP->entropy[p][rec][0]);
                nc_put_vara(ncid, KL_entropy_id, start, count, &resP->KL_entropy[p][rec][0]);
            }
            if (simP->mode_stats == 3)
            {
                for (int rec1=0; rec1<steps; rec1++)
                {
                    nc_put_vara(ncid, p_range_id, start, count, &resP->p_range[p][rec][rec1][0]);
                    nc_put_vara(ncid, p_exact_id, start, count, &resP->p_exact[p][rec][rec1][0]);
                    nc_put_vara(ncid, p_approx_id, start, count, &resP->p_approx[p][rec][rec1][0]);
                    nc_put_vara(ncid, cdf_theory_id, start, count, &resP->cdf_theory[p][rec][rec1][0]);
                }
            }
        }
    }

    nc_close(ncid);
    cout << "done!" << endl;
}

void write_prep_paras(int ncid, int *dimID, simulation *simP)
{
    int dim;
    for (unsigned i=0;i<simP->nVar;i++)
    {
        // cout << "preparing " << simP->vars[i].name << endl;
        nc_def_dim(ncid, strApp(simP->vars[i].name,"_dim").c_str(), simP->vars[i].steps, &dim);
        nc_def_var(ncid, simP->vars[i].name.c_str(), NC_DOUBLE, 1, &dim, dimID+i);
    }
}

void write_paras(int ncid, int *dimID, simulation *simP)
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
