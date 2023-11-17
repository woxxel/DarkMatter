#include <plog/Log.h>

void Simulation_Variable::initialize(double *modVarP, double *varP, unsigned varSz, string varName)
{
    /*
        modVarP (double*)
            pointer to variable value in model
        varP (double*)
            pointer to array containing iteration values of variable
        varSz (unsigned)
            number of iterations for this variable
        varName (string)
            Name of the variable

        TODO:
            * if additional variables are added, should check for consistency of 
                varSz, varName and array - or should it?
    */

    PLOG_DEBUG <<" Initializing " << varName << " with " << varSz << " variables";
    // cout << "current value of variable: " << *modVarP << endl;
    name = varName;

    iter = 0;           // value is the "iter"th-entry in provided array
    steps = varSz;      // total number of values provided in array

    valP = varP;        // store pointer to array

    paraP.push_back(modVarP);   // contains pointers to all iterable variables

    // // needs to be implemented as executed function when setting rate
    // if (mod.paras.drive == 1) // inhibitory rate determined by network parameters
    // {
    //     mod.paras.rate[0] = sim.rateWnt[sim.rateWnt_iter]*mod.paras.kappa*mod.paras.J_E[0]/mod.paras.J_I[0];
    //     // cout << " rate inhib: " << mod.paras.rate[0] << "\t rate exc: " << mod.paras.rate[1] << endl;
    // }

    // set initial values for all model parameters in paraP
    for (unsigned i=0;i<paraP.size();i++) {
        *paraP[i] = *valP;
    }
};

bool Simulation_Variable::iterate() // return false if end is reached, true else
{
    /*
        iterates variable by one (or resets if end is reached)
    */

    iter = (iter+1) % steps;
    for (unsigned i=0;i<paraP.size();i++) {
        *(paraP[i]) = *(valP + iter);
    }
    // cout << "Variable '" << name << "' value @ iter/steps: " << iter << "/" << steps << ": " << *(valP + iter) << endl;

    return iter>0;
};

void Simulation_Variable::print_status()
{
    PLOG_DEBUG << "status of variable " << name << ": \t" << "i=" << iter << "/" << steps << "\t val=" << *paraP[0];
};

void Simulation::print_simStatus()
{
    for (unsigned i=0;i<nVar;i++)
        vars[i].print_status();
}

void Simulation::reset_iteration() {
    PLOG_DEBUG << "resetting variables " << nVar;
    for (unsigned i=0;i<nVar;i++)
    {
        iter_status[i] = true;
        vars[i].iter = 0;
    }
};

void Simulation::set_vars(Model *modP, unsigned var_idx, unsigned l, unsigned p, unsigned s)
{
    if (order[var_idx].compare(0,3,"eta")==0) {
        vars[var_idx].initialize(&modP->layer[l].eta, &eta[0], eta.size(), "eta");
    } else if (order[var_idx].compare(0,3,"eps")==0) {
        vars[var_idx].initialize(&modP->layer[l].eps, &eps[0], eps.size(), "eps");
    } else if (order[var_idx].compare(0,7,"alpha_0")==0) {
        vars[var_idx].initialize(&modP->layer[l].population[p].alpha_0, &alpha_0[0], alpha_0.size(), "alpha_0");
    } else if (order[var_idx].compare(0,5,"Psi_0")==0) {
        vars[var_idx].initialize(&modP->layer[l].population[p].Psi_0, &Psi_0[0], Psi_0.size(), "Psi_0");
    } else if (order[var_idx].compare(0,7,"rateWnt")==0) {
        vars[var_idx].initialize(&modP->layer[l].population[p].rateWnt, &rateWnt[0], rateWnt.size(), "rateWnt");
    } else if (order[var_idx].compare(0,5,"tau_n")==0) {
        vars[var_idx].initialize(&modP->layer[l].population[p].tau_n, &tau_n[0], tau_n.size(), "tau_n");
    } else if (order[var_idx].compare(0,5,"tau_I")==0) {
        vars[var_idx].initialize(&modP->layer[l].population[p].psp[s].tau_I, &tau_I[0], tau_I.size(), "tau_I");
    } else if (order[var_idx].compare(0,7,"I_alpha")==0)
        vars[var_idx].initialize(&modP->I_alpha, &I_alpha[0], I_alpha.size(), "I_alpha");
    else if (order[var_idx].compare(0,6,"I_beta")==0)
        vars[var_idx].initialize(&modP->I_beta, &I_beta[0], I_beta.size(), "I_beta");
}

void Simulation::initialize(Model *modP)
{
    nVar = order.size();
    vars.resize(nVar);

    unsigned layer, population, psp;
    // cout << nVar << " variables in input " << endl;
    vector<int> l,p,s;

    for (unsigned i=0;i<nVar;i++)
    {
        // getting layers to change variables in from sim_pointer-input
        l.clear();
        l.push_back(sim_pointer[i][0]); // should be obtained from provided vector
        if (l[0]<0) {
            l.resize(modP->L);
            for (unsigned ll=0; ll<modP->L; ll++) l[ll] = ll;
        }
        for (unsigned ll=0; ll<l.size(); ll++) {
            layer = l[ll];

            // cout << "layer (" << ll << ")" << layer << endl;
            // getting populations to change variables in from sim_pointer-input
            p.clear();
            p.push_back(sim_pointer[i][1]);
            if (p[0]<0) {
                p.resize(modP->layer[layer].nPop);
                for (unsigned pp=0; pp<modP->layer[layer].nPop; pp++) p[pp] = pp;
            }
            for (unsigned pp=0; pp<p.size(); pp++) {
                population = p[pp];
                // cout << "population (" << pp << ")" << population << endl;

                // getting synapses to change variables in from sim_pointer-input
                s.clear();
                s.push_back(sim_pointer[i][2]);
                if (s[0]<0) {
                    s.resize(modP->layer[layer].population[population].nPSP);
                    for (unsigned ss=0; ss<modP->layer[layer].population[population].nPSP; ss++) s[ss] = ss;
                }
                for (unsigned ss=0; ss<s.size(); ss++) {
                    psp = s[ss];
                    // cout << "psp (" << ss << ") " << psp << endl;

                    // cout << "layer|population|synapse: " << layer << "|" << population << "|" << psp << endl;

                    set_vars(modP,i,layer,population,psp);
                }
            }
        }


        // cout << "this variable iterates layer|population|synapse: " << l << "|" << p << "|" << s << endl;


    }
    // cout << "init done "<< endl;
    steps = vars[0].steps;
    reset_iteration();
};

bool Simulation::run_iteration(Model *modP, Model *modP_approx)
{
    /*
        function manages the iteration through simulation variables

        variables are provided as input in the order in which they should be iterated: first one is the inner (faster) loop, last one the outer loop

        outer loops are only iterated, when inner variables have completed one loop (are marked as "iter_status=false" by "iterate()"-method)
    */
    unsigned loops=0;
    for (unsigned i=0;i<nVar;i++) {
        
        if (iter_status[i]) // only iterate variables tagged as "true" (to maintain hierarchy)
        {
            iter_status[i] = vars[i].iterate();
            vars[i].print_status();
            if (i>0)  // whenever an outer loop is iterated, all inner loops are reset to true
                iter_status[i-1] = true;
            
            // else if (vars[i].iter==0) { // reset axes when loop is at start
                // initiate_y_axis(modP);
                // initiate_y_axis(modP_approx);
            // }

            // if iteration is done and iter_status is still true, don't jump to next variable
            if (iter_status[i])
                break;
        }
        loops++;
    }
    // print_simStatus();
    
    return loops!=nVar;
};

// void Simulation::initiate_y_axis(Model *modP)
// {
//     Population_Simulation *popSimP;
//     // cout << "initiating y-axis..." << endl;
//     for (unsigned l = 0; l < modP->L; l++) {
//         for (unsigned p = 0; p < modP->layer[l].nPop; p++) {

//             popSimP = &modP->layer[l].population[p].simulation;

//             // popSimP->trans_DM_found = false;
//             // popSimP->trans_np_found = false;

//             // popSimP->trans_DM = NAN;
//             // popSimP->trans_np = NAN;
//         }
//     }

//     // modP->simulation.trans_inc = NAN;
//     // modP->simulation.trans_imp = NAN;

//     // modP->simulation.trans_inc_found = false;
//     // modP->simulation.trans_imp_found = false;
// }


void Simulation::store_results(Model * modP, Model * modP_approx)
{
    /*
        at the end of each iteration, write mode-appropriate values to results-file 
    */

//         unsigned a = resP->rate.size() - 1;
//         cout << "size = " << a << endl;
//         unsigned size = resP->rate[a].size()+1;
//         cout << "rate=" << paras.rate << " ,\t q=" << paras.q[0] << " ,\t alpha=" << paras.alpha[0] << " ,\t gamma=" << paras.gamma[0] << " ,\t chi=" << paras.chi[0] << endl;
        //! write: rate, q, alpha, alpha+alpha_0, sigma_V, gamma, chi, threshold transition, nu_no_peak, nu_inconsistent
    Population_Simulation *popSimP, *popSimP_approx;
    Population_Results *popResP, *popResP_approx;

    // cout << "rate: " << modP->layer[0].population[0].simulation.rate;
    // cout << ", alpha: " << modP->layer[0].population[0].alpha_0 << endl;
    for (unsigned l = 0; l < modP->L; l++) {
        for (unsigned p = 0; p < modP->layer[l].nPop; p++) {
            popSimP = &modP->layer[l].population[p].simulation;
            popResP = &modP->layer[l].population[p].results;

            popResP->rate[vars[1].iter][vars[0].iter] = popSimP->rate;
            popResP->q[vars[1].iter][vars[0].iter] = popSimP->q;

            popResP->gamma[vars[1].iter][vars[0].iter] = popSimP->gamma;
            popResP->delta[vars[1].iter][vars[0].iter] = popSimP->delta;
            popResP->rate_max[vars[1].iter][vars[0].iter] = popSimP->rate_max;

            popResP->chi[vars[1].iter][vars[0].iter] = popSimP->chi;
            popResP->I_balance[vars[1].iter][vars[0].iter] = popSimP->I_balance;

            if (popSimP->trans_DM_found) // && (!simP->trans_DM_found[l][p])
            {
                popResP->trans_DM[vars[1].iter].push_back(popSimP->trans_DM);
                modP->simulation.nTrans = max(popResP->trans_DM[vars[1].iter].size(),modP->simulation.nTrans);
            }

            if (popSimP->trans_np_found) // && (!simP->trans_np_found[l][p])
            {
                popResP->trans_np[vars[1].iter].push_back(popSimP->trans_np);
                modP->simulation.nTrans = max(popResP->trans_np[vars[1].iter].size(),modP->simulation.nTrans);
            }

            // writing approximate results
            if ((mode_stats == 2) || (mode_stats == 3))
            {
                popSimP_approx = &modP_approx->layer[l].population[p].simulation;
                popResP_approx = &modP_approx->layer[l].population[p].results;

                if (popSimP_approx->trans_DM_found)
                {
                    popResP_approx->trans_DM[vars[1].iter].push_back(popSimP_approx->trans_DM);
                    modP_approx->simulation.nTrans = max(popResP_approx->trans_DM[vars[1].iter].size(),modP_approx->simulation.nTrans);
                }

                if (popSimP_approx->trans_np_found)
                {
                    popResP_approx->trans_np[vars[1].iter].push_back(popSimP_approx->trans_np);
                    modP_approx->simulation.nTrans = max(popResP_approx->trans_np[vars[1].iter].size(),modP_approx->simulation.nTrans);
                }
            }

            if ((mode_stats == 0) || (mode_stats == 3))
            {
                popResP->regions[vars[1].iter][vars[0].iter] = popSimP->regions;
                if ((mode_stats == 2) || (mode_stats == 3))
                    popResP_approx->regions[vars[1].iter][vars[0].iter] = popSimP_approx->regions;
            }

            if (mode_stats == 1)
            {
                popResP->regions[vars[1].iter][vars[0].iter] = popSimP->regions;
                popResP->alpha_raw[vars[1].iter][vars[0].iter] = popSimP->alpha_raw;
                popResP->alpha[vars[1].iter][vars[0].iter] = popSimP->alpha;
                popResP->sigma_V[vars[1].iter][vars[0].iter] = popSimP->sigma_V;
            }

            if ((mode_stats == 2) || (mode_stats == 3))
            {
                popResP_approx->q[vars[1].iter][vars[0].iter] = popSimP_approx->q;
                popResP_approx->gamma[vars[1].iter][vars[0].iter] = popSimP_approx->gamma;
                popResP_approx->chi[vars[1].iter][vars[0].iter] = popSimP_approx->chi;

                popResP->KL_entropy[vars[1].iter][vars[0].iter] = popSimP->KL;
                popResP->entropy[vars[1].iter][vars[0].iter] = popSimP->entropy;
            }

            if (mode_stats == 4)
            {
                // for (unsigned z=0; z<infoParas.nZeta; z++)
                // {
                    // cout << "handing over data: " << modP->infoContent[p][z] << endl;
                popResP->infoContent[vars[1].iter][vars[0].iter] = popSimP->infoContent;
                // }
            }
        }
    }
    // cout << "inconsistent? " << modP->simulation.trans_inc_found << endl;
    // cout << "imp found: " << modP->simulation.trans_imp_found << ", " << modP->simulation.trans_imp << endl;


    if (modP->simulation.trans_imp_found) // && (!simP->trans_imp_found)
    {
        modP->results.trans_imp[vars[1].iter].push_back(modP->simulation.trans_imp);
        modP->simulation.nTrans = max(modP->results.trans_imp[vars[1].iter].size(),modP->simulation.nTrans);
        // cout << "imp transition found and storing" << endl;
    }
        // trans_imp_found = true;

    if (modP->simulation.trans_inc_found) // && (!simP->trans_inc_found)
    {
        modP->results.trans_inc[vars[1].iter].push_back(modP->simulation.trans_inc);
        modP->simulation.nTrans = max(modP->results.trans_inc[vars[1].iter].size(),modP->simulation.nTrans);
        // cout << "inc transition found and storing" << endl;
    }
        // trans_inc_found = true;

    if ((mode_stats == 2) || (mode_stats == 3)) {
        if (trans_imp_found_approx)
        {
            modP_approx->results.trans_imp[vars[1].iter].push_back(modP_approx->simulation.trans_imp);
            modP_approx->simulation.nTrans = max(modP_approx->results.trans_imp[vars[1].iter].size(),modP_approx->simulation.nTrans);
        }
        // trans_imp_found_approx = true;

        if (trans_inc_found_approx)
        {
            modP_approx->results.trans_inc[vars[1].iter].push_back(modP_approx->simulation.trans_inc);
            modP_approx->simulation.nTrans = max(modP_approx->results.trans_inc[vars[1].iter].size(),modP_approx->simulation.nTrans);
        }
        // trans_inc_found_approx = true;
    }
}
