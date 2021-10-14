#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

// #include "spdlog/spdlog.h"

#include "./computations/structures.h"
#include "./computations/functions.cpp"
#include "./computations/simFunctions.cpp"
#include "./io/readwrite.cpp"

int main(int argc, char** argv)
{

	// spdlog::set_level(spdlog::level::info);
	if ((argc < 4) or (argc > 4))
	{
		cout << "Please specify (only) the input and output file names!" << endl;
		return 1;
	}

//      parameters:
//          drive:
//          0 : all populations driven to have same firing rates
//          1 : only excitatory one driven
//          2 : recurrent, inhibitory population is driven by afferent spikes from excitatory population
	string fileModel = argv[1];
    string fileSim = argv[2];
	string fileOut = argv[3];

    simulation sim, *simP = &sim;
    model mod, *modP = &mod;
    results res, *resP = &res;

    read_model(fileModel,modP);
    read_simulation(fileSim,simP);

	sim.initialize(modP);

    initiate_results(modP,simP,resP);

	model mod_approx, *mod_approxP = &mod_approx;
	mod_approx = mod;

	cout << "calculate solutions for mode stats = " << sim.mode_stats << endl;

	while (sim.run_iteration(modP))
	{
		// cout << "start iteration... " << endl;

		// cout << "solving... " << endl;
		mod.solve_selfcon(sim.mode_calc);
		// cout << "done" << endl;
		// cout << "gamma : " << mod.paras.gamma[0] << endl;

		// cout << "solving approx... " << endl;
		if ((sim.mode_stats == 2) || (sim.mode_stats == 4))
		{
			mod_approx = mod;
			mod_approx.solve_selfcon(1);
		}
		// cout << "done" << endl;

		// cout << "comparing approx" << endl;
		if ((sim.mode_stats == 2) || (sim.mode_stats == 4))
			compare_approx(modP,mod_approxP,simP,resP);
		// cout << "done" << endl;

		// cout << "transitions... " << endl;
		find_transitions(modP,simP,resP);
		// cout << "done" << endl;

		// cout << "storing ..." << endl;
		sim.store_results(simP,modP,mod_approxP,resP);
		// cout << "done" << endl;
	};
	// cout << "computation done!" << endl;

	if (sim.steps == 1)
		write_theory(fileOut,resP);
	else {
        write_results(fileOut,simP,modP,resP);
    }
	return 0;
}
