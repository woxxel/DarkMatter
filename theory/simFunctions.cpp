#include <vector>

// #include "simFunctions.h"

void simulation_variable::initialize(double *varP, unsigned varSz, string varName)
{
    cout << "initializing var " << varName << endl;
    name = varName;
    iter = 0;
    max_iter = varSz;
    valP = varP;
};

bool simulation_variable::iterate() // return false if end is reached, true else
{
    val = *(valP + iter++);
    if (iter==max_iter)
    {
        iter=0;
        return false;
    } else {
        return true;
    }
};

void simulation_variable::print_status()
{
    cout << "status of variable " << name << ": \t" << "i=" << iter << "/" << max_iter << "\t val=" << val << endl;
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

void simulation::initialize()
{
    nVar = order.size();
    vars.resize(nVar);
    for (unsigned i=0;i<nVar;i++)
    {
        if (order[i].compare(0,7,"rateWnt")==0)
            vars[i].initialize(&rateWnt[0], rateWntSz, "rateWnt");
        else if (order[i].compare(0,5,"tau_G")==0)
            vars[i].initialize(&tau_G[0], tau_GSz, "tau_G");
        else if (order[i].compare(0,3,"eta")==0)
            vars[i].initialize(&eta[0], etaSz, "eta");
        else if (order[i].compare(0,3,"eps")==0)
            vars[i].initialize(&eps[0], epsSz, "eps");
        else if (order[i].compare(0,7,"alpha_0")==0)
            vars[i].initialize(&alpha_0[0], alpha_0Sz, "alpha_0");
        else if (order[i].compare(0,1,"n")==0)
            vars[i].initialize(&n[0], nSz, "n");
    }
    reset_iteration();
};

void simulation::run_iteration()
{
    unsigned loops;
    while (iterating) {
        loops = 0;
        for (unsigned i=0;i<nVar;i++) {
            if (iter_status[i])
            {
                iter_status[i] = vars[i].iterate();
                if (i>0)
                    iter_status[i-1] = true;
                if (iter_status[i])
                    break;
            }
            loops++;
        }
        print_simStatus();
        iterating = loops!=nVar;
    }
};
