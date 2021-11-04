
void simulation_variable::initialize(double *modVarP, unsigned modVarSz, double *varP, unsigned varSz, string varName)
{
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
