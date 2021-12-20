
void simulation_variable::initialize(double *modVarP, unsigned modVarSz, double *varP, unsigned varSz, string varName)
{
    cout <<" initializing " << varName << " with " << varSz << " variables" << endl;
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

    for (unsigned i=0;i<nVar;i++)
    {
        if (order[i].compare(0,7,"alpha_0")==0)
            vars[i].initialize(&modP->paras.alpha_0[0], modP->paras.Npop, &alpha_0[0], alpha_0Sz, "alpha_0");
        else if (order[i].compare(0,7,"rateWnt")==0)
            vars[i].initialize(&modP->paras.rate[0], modP->paras.Npop, &rateWnt[0], rateWntSz, "rateWnt");
        else if (order[i].compare(0,5,"tau_G")==0)
            vars[i].initialize(&modP->paras.tau_G, 1, &tau_G[0], tau_GSz, "tau_G");
        else if (order[i].compare(0,3,"eta")==0)
            vars[i].initialize(&modP->paras.eta, 1, &eta[0], etaSz, "eta");
        else if (order[i].compare(0,3,"eps")==0)
            vars[i].initialize(&modP->paras.eps, 1, &eps[0], epsSz, "eps");
        else if (order[i].compare(0,1,"n")==0)
            vars[i].initialize(&modP->paras.n, 1, &n[0], nSz, "n");
        else if (order[i].compare(0,4,"zeta")==0)
            vars[i].initialize(&modP->paras.zeta, 1, &zeta[0], zetaSz, "zeta");
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
    // print_simStatus();
    modP->set_weights(); // could rather be implemented via setter to according variables
    return loops!=nVar;
};

void simulation::initiate_y_axis(model *modP)
{
    // cout << "initiating y-axis..." << endl;
    for (unsigned p = 0; p < modP->paras.Npop; p++)
    {
        trans_DM_found[p] = false;
        trans_np_found[p] = false;

        modP->trans_DM[p] = NAN;
        modP->trans_np[p] = NAN;
        modP->trans_inc[p] = NAN;
        modP->trans_imp[p] = NAN;
        trans_DM_found_approx[p] = false;
        trans_np_found_approx[p] = false;
    }

    trans_inc_found = false;
    trans_imp_found = false;

    trans_inc_found_approx = false;
    trans_imp_found_approx = false;
}


void simulation::store_results(model * modP, model * mod_approxP, results * resP)
{
//         unsigned a = resP->rate.size() - 1;
//         cout << "size = " << a << endl;
//         unsigned size = resP->rate[a].size()+1;
//         cout << "rate=" << paras.rate << " ,\t q=" << paras.q[0] << " ,\t alpha=" << paras.alpha[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
        //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent
    for (unsigned p = 0; p < modP->paras.Npop; p++)
    {
        resP->rate[p][vars[1].iter][vars[0].iter] = modP->paras.rate[p];
        resP->q[p][vars[1].iter][vars[0].iter] = modP->paras.q[p];

        resP->gamma[p][vars[1].iter][vars[0].iter] = modP->paras.gamma[p];
        resP->chi[p][vars[1].iter][vars[0].iter] = modP->paras.chi[p];
        resP->delta[p][vars[1].iter][vars[0].iter] = modP->paras.delta[p];
        resP->I_balance[p][vars[1].iter][vars[0].iter] = modP->paras.I_balance[p];

        if ((!trans_DM_found[p]) && (modP->trans_DM_found[p]))
        {
            resP->trans_DM[p][vars[1].iter] = modP->trans_DM[p];
            trans_DM_found[p] = true;
        }
        if ((!trans_np_found[p]) && (modP->trans_np_found[p]))
        {
            resP->trans_np[p][vars[1].iter] = modP->trans_np[p];
            trans_np_found[p] = true;
        }
        if ((!trans_imp_found) && (modP->trans_imp_found))
        {
            resP->trans_imp[p][vars[1].iter] = modP->trans_imp[p];
            trans_imp_found = true;
        }
        if ((!trans_inc_found) && (modP->trans_inc_found))
        {
            resP->trans_inc[p][vars[1].iter] = modP->trans_inc[p];
            trans_inc_found = true;
        }
        if ((mode_stats == 2) || (mode_stats == 3))
        {
            if ((!trans_DM_found_approx[p]) && (mod_approxP->trans_DM_found[p]))
            {
                resP->trans_DM_approx[p][vars[1].iter] = mod_approxP->trans_DM[p];
                trans_DM_found_approx[p] = true;
            }
            if ((!trans_np_found_approx[p]) && (mod_approxP->trans_np_found[p]))
            {
                resP->trans_np_approx[p][vars[1].iter] = mod_approxP->trans_np[p];
                trans_np_found_approx[p] = true;
            }
            if ((!trans_imp_found_approx) && (mod_approxP->trans_imp_found))
            {
                resP->trans_imp_approx[p][vars[1].iter] = mod_approxP->trans_imp[p];
                trans_imp_found_approx = true;
            }
            if ((!trans_inc_found_approx) && (mod_approxP->trans_inc_found))
            {
                resP->trans_inc_approx[p][vars[1].iter] = mod_approxP->trans_inc[p];
                trans_inc_found_approx = true;
            }
        }

        if ((mode_stats == 0) || (mode_stats == 3))
        {
            resP->regions[p][vars[1].iter][vars[0].iter] = modP->paras.regions[p];
            if ((mode_stats == 2) || (mode_stats == 3))
                resP->regions_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.regions[p];
        }

        if (mode_stats == 1)
        {
            resP->alpha_raw[p][vars[1].iter][vars[0].iter] = modP->paras.alpha_raw[p];
            resP->alpha[p][vars[1].iter][vars[0].iter] = modP->paras.alpha[p];
            resP->sigma_V[p][vars[1].iter][vars[0].iter] = modP->paras.sigma_V[p];
        }

        if ((mode_stats == 2) || (mode_stats == 3))
        {
            resP->q_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.q[p];
            resP->gamma_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.gamma[p];
            resP->chi_approx[p][vars[1].iter][vars[0].iter] = mod_approxP->paras.chi[p];

            resP->KL_entropy[p][vars[1].iter][vars[0].iter] = modP->paras.KL[p];
            resP->entropy[p][vars[1].iter][vars[0].iter] = modP->paras.entropy[p];
        }

        if (mode_stats == 4)
        {
            // for (unsigned z=0; z<infoParas.nZeta; z++)
            // {
                // cout << "handing over data: " << modP->infoContent[p][z] << endl;
            resP->infoContent[p][vars[1].iter][vars[0].iter] = modP->infoContent[p];
            // }
        }
    }
}
