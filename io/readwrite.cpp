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

void write_to_ncid(int ncid, const char *varName, nc_type varType, int dimSz, const int *dimids, double *varData, hyperslab hslab)
{
    int varid;
    nc_def_var(ncid, varName, varType, dimSz, dimids, &varid);
    cout << "data to write: " << *(varName+0) << *(varName+1) << *(varName+2) << *(varName+3) << endl;

    double var;
    for (int i=0; i<10; i++) {
        var = *(varData + i);
        cout << "data (" << i << "): " << var << endl;
    }
    // nc_enddef(ncid);

    nc_put_vars(ncid, varid, &hslab.startp[0], &hslab.countp[0], &hslab.stridep[0], varData);
};

void read_model(string fileModel, model *modP)
{
    cout << "reading model parameters from " << fileModel << "... ";
    int ncid;

    nc_open(fileModel.c_str(), NC_NOWRITE, &ncid);

    get_from_ncid(ncid, "Npop", &modP->paras.Npop);
    modP->resize();

    // get the membrane constants
    get_from_ncid(ncid, "tau_A", &modP->paras.tau_A);
    get_from_ncid(ncid, "tau_N", &modP->paras.tau_N);
    get_from_ncid(ncid, "tau_M", &modP->paras.tau_M);
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

    // read variables from ncid
    get_from_ncid(ncid, "mode_calc", &simP->mode_calc);
    get_from_ncid(ncid, "mode_stats", &simP->mode_stats);
    // cout << "mode calc: " << simP->mode_calc << endl;
    // cout << "mode stats: " << simP->mode_stats << endl;
    // cout << "sizes (in Sim) - n: " << simP->nSz << ", alpha_0: " << simP->alpha_0Sz << ", tau_G: " << simP->tau_GSz << ", rate: " << simP->rateWntSz << endl;
    get_from_ncid(ncid, "n", &simP->n[0]);
    get_from_ncid(ncid, "alpha_0", &simP->alpha_0[0]);
    get_from_ncid(ncid, "tau_G", &simP->tau_G[0]);
    get_from_ncid(ncid, "rateWnt", &simP->rateWnt[0]);
    get_from_ncid(ncid, "eps", &simP->eps[0]);
    get_from_ncid(ncid, "eta", &simP->eta[0]);
    
    if (simP->nSz == 1)
        simP->n.resize(simP->steps,simP->n[0]);

    if (simP->alpha_0Sz == 1)
        simP->alpha_0.resize(simP->steps,simP->alpha_0[0]);

    if (simP->tau_GSz == 1)
        simP->tau_G.resize(simP->steps,simP->tau_G[0]);

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


void write_sharks(string fileOut, simulation sim, model mod, results *resP)
{
	cout << "writing shark data to file '" << fileOut << "'... ";

    int ncid, dimids[3], steps_dim, Npop_dim;
    int steps = resP->steps;
    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    nc_def_dim(ncid, "Npop", mod.paras.Npop, &Npop_dim);
    nc_def_dim(ncid, "steps", steps, &steps_dim);

    dimids[0] = Npop_dim;
    dimids[1] = steps_dim;
    dimids[2] = steps_dim;

    int DM_id, np_id, inc_id, imp_id;
    int gamma_id, chi_id, regions_id;
    int gamma_approx_id, chi_approx_id, entropy_id, KL_entropy_id;
    nc_def_var(ncid, "gamma", NC_DOUBLE, 3, &dimids[0], &gamma_id);
    nc_def_var(ncid, "chi", NC_DOUBLE, 3, &dimids[0], &chi_id);
    nc_def_var(ncid, "regions", NC_DOUBLE, 3, &dimids[0], &regions_id);

    if (sim.mode_stats == 4)
    {
        nc_def_var(ncid, "gamma_approx", NC_DOUBLE, 3, &dimids[0], &gamma_approx_id);
        nc_def_var(ncid, "chi_approx", NC_DOUBLE, 3, &dimids[0], &chi_approx_id);
        nc_def_var(ncid, "entropy", NC_DOUBLE, 3, &dimids[0], &entropy_id);
        nc_def_var(ncid, "KL_entropy", NC_DOUBLE, 3, &dimids[0], &KL_entropy_id);
    }


    nc_def_var(ncid, "DM_trans", NC_DOUBLE, 2, &dimids[0], &DM_id);
    nc_def_var(ncid, "np_trans", NC_DOUBLE, 2, &dimids[0], &np_id);
    nc_def_var(ncid, "inc_trans", NC_DOUBLE, 2, &dimids[0], &inc_id);
    nc_def_var(ncid, "imp_trans", NC_DOUBLE, 2, &dimids[0], &imp_id);

    nc_enddef(ncid);

    size_t start[3], count[3];
    start[2] = 0;
    count[0] = 1;
    count[2] = steps;

    for (int p=0; p<mod.paras.Npop; p++)
    {
        start[0] = p;

        start[1] = 0;
        count[1] = steps;
        nc_put_vara(ncid, DM_id, start, count, &resP->trans_DM[p][0]);
        nc_put_vara(ncid, np_id, start, count, &resP->trans_np[p][0]);
        nc_put_vara(ncid, inc_id, start, count, &resP->trans_inc[p][0]);
        nc_put_vara(ncid, imp_id, start, count, &resP->trans_imp[p][0]);

        count[1] = 1;
        for (int rec=0; rec<steps; rec++)
        {
            start[1] = rec;
            nc_put_vara(ncid, gamma_id, start, count, &resP->gamma[p][rec][0]);
            nc_put_vara(ncid, chi_id, start, count, &resP->chi[p][rec][0]);
            nc_put_vara(ncid, regions_id, start, count, &resP->regions[p][rec][0]);

            if (sim.mode_stats == 4)
            {
                nc_put_vara(ncid, gamma_approx_id, start, count, &resP->gamma_approx[p][rec][0]);
                nc_put_vara(ncid, chi_approx_id, start, count, &resP->chi_approx[p][rec][0]);
                nc_put_vara(ncid, entropy_id, start, count, &resP->entropy[p][rec][0]);
                nc_put_vara(ncid, KL_entropy_id, start, count, &resP->KL_entropy[p][rec][0]);
            }
        }
    }

    nc_close(ncid);
}

void write_stats(string fileOut, model mod, simulation sim, results * resP)
{
    cout << "writing simulation data to file '" << fileOut << "'..." << endl;

    int ncid, dimids[2], Nc_dim1, Nc_dim2, alpha_0_dim;
    nc_create(fileOut.c_str(), NC_CLOBBER, &ncid);

    int dim1 = sim.max_ax[1];
    int dim2 = sim.max_ax[0];

    nc_def_dim(ncid, "dim1", dim1, &Nc_dim1);
    nc_def_dim(ncid, "dim2", dim2, &Nc_dim2);

    dimids[0] = Nc_dim1;
    dimids[1] = Nc_dim2;

    write_to_ncid(ncid,"gamma", NC_DOUBLE, 2, &dimids[0], &resP->gamma[0][0][0]);


    if (sim.mode_stats == 1)
    {
        if (sim.eps.size() == 1)
            sim.eps.resize(dim2,sim.eps[0]);
        if (sim.eta.size() == 1)
            sim.eta.resize(dim2,sim.eta[0]);
        if (sim.n.size() == 1)
            sim.n.resize(dim2,sim.n[0]);
        if (sim.tau_G.size() == 1)
            sim.tau_G.resize(dim2,sim.tau_G[0]);

//                 write_stats_mode_1(&writeResults,resP,0,"",sim.max_ax[1],Nc_dim1,sim.max_ax[0],Nc_dim2);
//                 if sim.eps.size() > 1:
        cout << "now write" << endl;

        cout << "eps size: " << sim.eps.size() << endl;
        cout << "eta size: " << sim.eta.size() << endl;
        cout << "n size: " << sim.n.size() << endl;
        cout << "tau_G size: " << sim.tau_G.size() << endl;

        write_to_ncid(ncid,"eps", NC_DOUBLE, 2, &dimids[0], &sim.eps[0]);
        write_to_ncid(ncid,"eta", NC_DOUBLE, 2, &dimids[0], &sim.eta[0]);
        write_to_ncid(ncid,"n", NC_DOUBLE, 2, &dimids[0], &sim.n[0]);
        write_to_ncid(ncid,"tau_G", NC_DOUBLE, 2, &dimids[0], &sim.tau_G[0]);


        // NcVar *eps = writeResults.add_var("eps", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *eta = writeResults.add_var("eta", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *n = writeResults.add_var("n", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *tau_G = writeResults.add_var("tau_G", ncDouble, Nc_dim1, Nc_dim2);

        hyperslab hslab_1 = { {0,0,0}, {1,size_t(dim1),size_t(dim2)}, {1,1,1} };
        write_to_ncid(ncid,"rateWnt", NC_DOUBLE, 2, &dimids[0], &resP->rate[0][0][0], hslab_1);
        write_to_ncid(ncid,"q", NC_DOUBLE, 2, &dimids[0], &resP->q[0][0][0], hslab_1);
        write_to_ncid(ncid,"alpha_raw", NC_DOUBLE, 2, &dimids[0], &resP->alpha_raw[0][0][0], hslab_1);
        write_to_ncid(ncid,"alpha", NC_DOUBLE, 2, &dimids[0], &resP->alpha[0][0][0], hslab_1);
        write_to_ncid(ncid,"sigma_V", NC_DOUBLE, 2, &dimids[0], &resP->sigma_V[0][0][0], hslab_1);
        write_to_ncid(ncid,"gamma", NC_DOUBLE, 2, &dimids[0], &resP->gamma[0][0][0], hslab_1);
        write_to_ncid(ncid,"I_balance", NC_DOUBLE, 2, &dimids[0], &resP->I_balance[0][0][0], hslab_1);
        write_to_ncid(ncid,"chi", NC_DOUBLE, 2, &dimids[0], &resP->chi[0][0][0], hslab_1);

        // NcVar *rateWnt = writeResults.add_var("rateWnt", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *q = writeResults.add_var("q", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *alpha_raw = writeResults.add_var("alpha_raw", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *alpha = writeResults.add_var("alpha", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *sigma_V = writeResults.add_var("sigma_V", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *gamma = writeResults.add_var("gamma", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *I_balance = writeResults.add_var("I_balance", ncDouble, Nc_dim1, Nc_dim2);
        // NcVar *chi = writeResults.add_var("chi", ncDouble, Nc_dim1, Nc_dim2);


        // for (int rec=0; rec<dim1; ++rec)
        // {
        //
        //     // eps->put_rec(&sim.eps[0],rec);
        //     // eta->put_rec(&sim.eta[0],rec);
        //     // n->put_rec(&sim.n[0],rec);
        //     // tau_G->put_rec(&sim.tau_G[0],rec);
        //
        //     rateWnt->put_rec(&resP->rate[0][rec][0],rec);
        //     q->put_rec(&resP->q[0][rec][0],rec);
        //     alpha_raw->put_rec(&resP->alpha_raw[0][rec][0],rec);
        //     alpha->put_rec(&resP->alpha[0][rec][0],rec);
        //     sigma_V->put_rec(&resP->sigma_V[0][rec][0],rec);
        //     gamma->put_rec(&resP->gamma[0][rec][0],rec);
        //     I_balance->put_rec(&resP->I_balance[0][rec][0],rec);
        //     chi->put_rec(&resP->chi[0][rec][0],rec);
        // }
        write_to_ncid(ncid,"nu_imp", NC_DOUBLE, 1, &Nc_dim1, &resP->trans_imp[0][0]);

        // NcVar *trans_imp = writeResults.add_var("nu_implausible", ncDouble, Nc_dim1);
        // trans_imp->put(&resP->trans_imp[0][0], dim1);

        if (mod.paras.Npop == 2)
        {
            cout << "#2" << endl;
            hyperslab hslab_2 = { {1,0,0}, {1,size_t(dim1),size_t(dim2)}, {1,1,1} };

            write_to_ncid(ncid,"rateWnt_exc", NC_DOUBLE, 2, &dimids[0], &resP->rate[1][0][0], hslab_2);
            write_to_ncid(ncid,"q_exc", NC_DOUBLE, 2, &dimids[0], &resP->q[1][0][0], hslab_2);
            write_to_ncid(ncid,"alpha_raw_exc", NC_DOUBLE, 2, &dimids[0], &resP->alpha_raw[1][0][0], hslab_1);
            write_to_ncid(ncid,"alpha_exc", NC_DOUBLE, 2, &dimids[0], &resP->alpha[1][0][0], hslab_2);
            write_to_ncid(ncid,"sigma_V_exc", NC_DOUBLE, 2, &dimids[0], &resP->sigma_V[1][0][0], hslab_2);
            write_to_ncid(ncid,"gamma_exc", NC_DOUBLE, 2, &dimids[0], &resP->gamma[1][0][0], hslab_2);
            write_to_ncid(ncid,"I_balance_exc", NC_DOUBLE, 2, &dimids[0], &resP->I_balance[1][0][0], hslab_2);
            write_to_ncid(ncid,"chi_exc", NC_DOUBLE, 2, &dimids[0], &resP->chi[1][0][0], hslab_2);

            // NcVar *rateWnt = writeResults.add_var("rateWnt_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *q = writeResults.add_var("q_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *alpha_raw = writeResults.add_var("alpha_raw_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *alpha = writeResults.add_var("alpha_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *sigma_V = writeResults.add_var("sigma_V_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *gamma = writeResults.add_var("gamma_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *I_balance = writeResults.add_var("I_balance_exc", ncDouble, Nc_dim1, Nc_dim2);
            // NcVar *chi = writeResults.add_var("chi_exc", ncDouble, Nc_dim1, Nc_dim2);
            // cout << "now write" << endl;
            //
            // for (int rec=0; rec<dim1; ++rec)
            // {
            //     rateWnt->put_rec(&resP->rate[1][rec][0],rec);
            //     q->put_rec(&resP->q[1][rec][0],rec);
            //     alpha_raw->put_rec(&resP->alpha_raw[1][rec][0],rec);
            //     alpha->put_rec(&resP->alpha[1][rec][0],rec);
            //     sigma_V->put_rec(&resP->sigma_V[1][rec][0],rec);
            //     gamma->put_rec(&resP->gamma[1][rec][0],rec);
            //     I_balance->put_rec(&resP->I_balance[1][rec][0],rec);
            //     chi->put_rec(&resP->chi[1][rec][0],rec);
            // }

            write_to_ncid(ncid,"nu_imp_exc", NC_DOUBLE, 1, &Nc_dim1, &resP->trans_imp[1][0]);
            // NcVar *trans_imp = writeResults.add_var("nu_implausible_exc", ncDouble, Nc_dim1);
            // trans_imp->put(&resP->trans_imp[1][0], dim1);


//                     write_stats_mode_1(&writeResults,resP,1,"_exc",sim.max_ax[1],Nc_dim1,sim.max_ax[0],Nc_dim2);
        }
    }

    if (sim.mode_stats == 2)
    {
        nc_def_dim(ncid, "alpha_0Sz", sim.alpha_0Sz, &alpha_0_dim);
        // write_to_ncid(ncid,"alpha_0Sz", NC_DOUBLE, 1, &Nc_dim1, &resP->trans_imp[0]);
        // write_stats_mode_2(&writeResults,alpha_0_dim,resP,0,"");
        write_stats_mode_2(ncid,alpha_0_dim,resP,0,"");

        if (mod.paras.Npop == 2)
            // write_stats_mode_2(&writeResults,alpha_0_dim,resP,1,"_exc");
            write_stats_mode_2(ncid,alpha_0_dim,resP,1,"_exc");
    }
    nc_close(ncid);
}
//         if (sim.mode_stats == 3)
//         {
// //                 int p_Sz = resP->steps;ctor<double> regions;
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



void write_stats_mode_2(int ncid, const int dimid, results *resP, int p, string addStr)
{
//         cout << "#2" << endl;
//     cout << strApp("inc_trans",addStr) << endl;
    write_to_ncid(ncid,"DM_trans", NC_DOUBLE, 1, &dimid, &resP->trans_DM[p][0]);
    write_to_ncid(ncid,"no_peak_trans", NC_DOUBLE, 1, &dimid, &resP->trans_np[p][0]);
    write_to_ncid(ncid,"inc_trans", NC_DOUBLE, 1, &dimid, &resP->trans_inc[p][0]);
    write_to_ncid(ncid,"implausible_trans", NC_DOUBLE, 1, &dimid, &resP->trans_imp[p][0]);
    // nc_close(ncid);

//         NcVar *trans_DM = writeResults->add_var("DM_trans", ncDouble, Nc_dim);
//         NcVar *trans_np = writeResults->add_var("no_peak_trans", ncDouble, Nc_dim);
//         NcVar *trans_inc = writeResults->add_var("inc_trans", ncDouble, Nc_dim);
//         NcVar *trans_imp = writeResults->add_var("nu_implausible", ncDouble, Nc_dim);
//
// //         cout << resP->trans_DM.size() << endl;
// //         cout << resP->trans_inc.size() << endl;
//         trans_DM->put(&resP->trans_DM[p][0], resP->steps);
//         trans_np->put(&resP->trans_np[p][0], resP->steps);
//         trans_inc->put(&resP->trans_inc[p][0], resP->steps);
//         trans_imp->put(&resP->trans_imp[p][0], resP->steps);
}

string strApp(string str1, string str2)
{
        string str="";
        str.append(str1);
        str.append(str2);
        return str;
}
