
void simulation_variable::initialize(double *modVarP, unsigned modVarSz, double *varP, unsigned varSz, string varName)
{
    /*
        modVarP (double*)
            pointer to variable value in model
        modVarSz (unsigned)
            number of different variables within the network (=Npop)
        varP (double*)
            pointer to array containing iteration values of variable
        varSz (unsigned)
            number of iterations for this variable
        varName (string)
            Name of the variable
    */
    cout <<" initializing " << varName << " with " << varSz << " variables" << endl;
    cout << "current value of variable: " << *modVarP << endl;
    name = varName;

    iter = 0;
    steps = varSz;

    valP = varP;

    paraP = modVarP;
    paraSz = modVarSz;

    // // needs to be implemented as executed function when setting rate
    // if (mod.paras.drive == 1) // inhibitory rate determined by network parameters
    // {
    //     mod.paras.rate[0] = sim.rateWnt[sim.rateWnt_iter]*mod.paras.kappa*mod.paras.J_E[0]/mod.paras.J_I[0];
    //     // cout << " rate inhib: " << mod.paras.rate[0] << "\t rate exc: " << mod.paras.rate[1] << endl;
    // }

    for (unsigned i=0;i<paraSz;i++)
        *(paraP+i) = *valP;
};

bool simulation_variable::iterate() // return false if end is reached, true else
{
    iter = (iter+1) % steps;
    // cout << "here name:" << name << ", iter/steps: " << iter << "/" << steps << endl;
    // cout << *(valP + iter) << endl;
    for (unsigned i=0;i<paraSz;i++)
        *(paraP+i) = *(valP + iter);

    // *onUpdate();

    return iter>0;
};

void simulation_variable::print_status()
{
    cout << "status of variable " << name << ": \t" << "i=" << iter << "/" << steps << "\t val=" << *paraP << endl;
};

void simulation::print_simStatus()
{
    for (unsigned i=0;i<nVar;i++)
    {
        cout << vars[i].name << endl;
        vars[i].print_status();
    }
}

void simulation::reset_iteration() {
    for (unsigned i=0;i<nVar;i++)
    {
        iter_status[i] = true;
        vars[i].iter = 0;
    }
};

void simulation::initialize(model *modP)
{
    nVar = order.size();
    vars.resize(nVar);

    cout << nVar << " variables in input " << endl;

    for (unsigned i=0;i<nVar;i++)
    {
        int l,p,s;
        l=sim_pointer[i][0]; // should be obtained from provided vector
        p=sim_pointer[i][1];
        s=sim_pointer[i][2];

        cout << "this variable iterates layer|population|synapse: " << l << "|" << p << "|" << s << endl;

        if (order[i].compare(0,3,"eta")==0)
            vars[i].initialize(&modP->layer[l].eta, 1, &eta[0], eta.size(), "eta");
        else if (order[i].compare(0,3,"eps")==0)
            vars[i].initialize(&modP->layer[l].eps, 1, &eps[0], eps.size(), "eps");
        else if (order[i].compare(0,7,"alpha_0")==0)
            vars[i].initialize(&modP->layer[l].population[p].alpha_0, 1, &alpha_0[0], alpha_0.size(), "alpha_0");
        else if (order[i].compare(0,7,"rateWnt")==0)
            vars[i].initialize(&modP->layer[l].population[p].rateWnt, 1, &rateWnt[0], rateWnt.size(), "rateWnt");
        else if (order[i].compare(0,5,"tau_n")==0)
            vars[i].initialize(&modP->layer[l].population[p].tau_n, 1, &tau_n[0], tau_n.size(), "tau_n");
        else if (order[i].compare(0,5,"tau_I")==0)
            vars[i].initialize(&modP->layer[l].population[p].psp[s].tau_I, 1, &tau_I[0], tau_I.size(), "tau_I");
        else if (order[i].compare(0,7,"I_alpha")==0)
            vars[i].initialize(&modP->paras.I_alpha, 1, &I_alpha[0], I_alpha.size(), "I_alpha");
        else if (order[i].compare(0,6,"I_beta")==0)
            vars[i].initialize(&modP->paras.I_beta, 1, &I_beta[0], I_beta.size(), "I_beta");
    }
    steps = vars[0].steps;
    reset_iteration();
};

bool simulation::run_iteration(model *modP)
{
    unsigned loops=0;
    for (unsigned i=0;i<nVar;i++) {
        if (iter_status[i])
        {
            // cout << "iterate: " << i << endl;
            iter_status[i] = vars[i].iterate();
            if (i>0)
                iter_status[i-1] = true;
            else if (vars[i].iter==0)
            {
                initiate_y_axis(modP);
            }
            if (iter_status[i])
                break;
        }
        loops++;
    }
    print_simStatus();
    modP->set_weights(); // could rather be implemented via setter to according variables
    return loops!=nVar;
};

void simulation::initiate_y_axis(model *modP)
{
    // cout << "initiating y-axis..." << endl;
    for (unsigned l = 0; l < modP->L; l++) {

        for (unsigned p = 0; p < modP->paras.Npop; p++)
        {
            modP->layer[l].population[p].simulation.trans_DM_found = false;
            modP->layer[l].population[p].simulation.trans_np_found = false;

            // modP->layer[l].population[p].trans_DM[p] = NAN;
            // modP->layer[l].population[p].trans_np[p] = NAN;
            // modP->layer[l].population[p].trans_inc[p] = NAN;
            // modP->layer[l].population[p].trans_imp[p] = NAN;
            modP->layer[l].population[p].simulation.trans_DM_found_approx = false;
            modP->layer[l].population[p].simulation.trans_np_found_approx = false;
        }
    }

    trans_inc_found = false;
    trans_imp_found = false;

    trans_inc_found_approx = false;
    trans_imp_found_approx = false;
}


void simulation::store_results(model * modP, model * mod_approxP)
{
// //         unsigned a = resP->rate.size() - 1;
// //         cout << "size = " << a << endl;
// //         unsigned size = resP->rate[a].size()+1;
// //         cout << "rate=" << paras.rate << " ,\t q=" << paras.q[0] << " ,\t alpha=" << paras.alpha[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
//         //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent
//
//
//
//     for (unsigned p = 0; p < modP->paras.Npop; p++)
//     {
//         resP->rate[p][vars[1].iter][vars[0].iter] = modP->paras.rate[p];
//         resP->q[p][vars[1].iter][vars[0].iter] = modP->paras.q[p];
//
//         resP->gamma[p][vars[1].iter][vars[0].iter] = modP->paras.gamma[p];
//         resP->chi[p][vars[1].iter][vars[0].iter] = modP->paras.chi[p];
//         resP->delta[p][vars[1].iter][vars[0].iter] = modP->paras.delta[p];
//         resP->I_balance[p][vars[1].iter][vars[0].iter] = modP->paras.I_balance[p];
//
//         if ((!trans_DM_found[p]) && (modP->trans_DM_found[p]))
//         {
//             resP->trans_DM[p][vars[1].iter] = modP->trans_DM[p];
//             trans_DM_found[p] = true;
//         }
//         if ((!trans_np_found[p]) && (modP->trans_np_found[p]))
//         {
//             resP->trans_np[p][vars[1].iter] = modP->trans_np[p];
//             trans_np_found[p] = true;
//         }
//         if ((!trans_imp_found) && (modP->trans_imp_found))
//         {
//             resP->trans_imp[p][vars[1].iter] = modP->trans_imp[p];
//             trans_imp_found = true;
//         }
//         if ((!trans_inc_found) && (modP->trans_inc_found))
//         {
//             resP->trans_inc[p][vars[1].iter] = modP->trans_inc[p];
//             trans_inc_found = true;
//         }
//         if ((mode_stats == 2) || (mode_stats == 3))
//         {
//             if ((!trans_DM_found_approx[p]) && (mod_approxP->trans_DM_found[p]))
//             {
//                 resP->trans_DM_approx[p][vars[1].iter] = mod_approxP->trans_DM[p];
//                 trans_DM_found_approx[p] = true;
//             }
//             if ((!trans_np_found_approx[p]) && (mod_approxP->trans_np_found[p]))
//             {
//                 resP->trans_np_approx[p][vars[1].iter] = mod_approxP->trans_np[p];
//                 trans_np_found_approx[p] = true;
//             }
//             if ((!trans_imp_found_approx) && (mod_approxP->trans_imp_found))
//             {
//                 resP->trans_imp_approx[p][vars[1].iter] = mod_approxP->trans_imp[p];
//                 trans_imp_found_approx = true;
//             }
//             if ((!trans_inc_found_approx) && (mod_approxP->trans_inc_found))
//             {
//                 resP->trans_inc_approx[p][vars[1].iter] = mod_approxP->trans_inc[p];
//                 trans_inc_found_approx = true;
//             }
//         }
//
//         if ((mode_stats == 0) || (mode_stats == 3))
//         {
//             resP->regions[p][vars[1].iter][vars[0].iter] = modP->paras.regions[p];
//             if ((mode_stats == 2) || (mode_stats == 3))
//                 resP->regions_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.regions[p];
//         }
//
//         if (mode_stats == 1)
//         {
//             resP->alpha_raw[p][vars[1].iter][vars[0].iter] = modP->paras.alpha_raw[p];
//             resP->alpha[p][vars[1].iter][vars[0].iter] = modP->paras.alpha[p];
//             resP->sigma_V[p][vars[1].iter][vars[0].iter] = modP->paras.sigma_V[p];
//         }
//
//         if ((mode_stats == 2) || (mode_stats == 3))
//         {
//             resP->q_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.q[p];
//             resP->gamma_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.gamma[p];
//             resP->chi_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.chi[p];
//
//             resP->KL_entropy[p][vars[1].iter][vars[0].iter] = modP->paras.KL[p];
//             resP->entropy[p][vars[1].iter][vars[0].iter] = modP->paras.entropy[p];
//         }
//
//         if (mode_stats == 4)
//         {
//             // for (unsigned z=0; z<infoParas.nZeta; z++)
//             // {
//                 // cout << "handing over data: " << modP->infoContent[p][z] << endl;
//             resP->infoContent[p][vars[1].iter][vars[0].iter] = modP->infoContent[p];
//             // }
//         }
//     }
}
