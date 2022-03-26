// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
// #include "spdlog/spdlog.h"

#include "./computations/structures.h"
#include "./computations/functions.h"
#include "./computations/functions.cpp"
#include "./computations/modFunctions.cpp"
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

	// parameters:
	//  	drive:
	//  	0 : all populations driven to have same firing rates
	//  	1 : only excitatory one driven
	//  	2 : recurrent, inhibitory population is driven by afferent spikes from excitatory population
	string fileModel = argv[1];
    string fileSim = argv[2];
	string fileOut = argv[3];

    Simulation sim, *simP = &sim;
    Model mod, *modP = &mod;
	Model mod_approx, *modP_approx = &mod_approx;
    // Results res, *resP = &res;

    read_model(fileModel,modP);
    read_model(fileModel,modP_approx);

    read_simulation(fileSim,simP);

	sim.initialize(modP);
	sim.initialize(modP_approx);

    initiate_results(modP,simP);
    initiate_results(modP_approx,simP);

	// Model *modP_approx = new Model(modP);
	// mod_approx = mod;

	cout << "calculate solutions for mode stats = " << sim.mode_stats << endl;

	while (sim.run_iteration(modP,modP_approx))
	{
		// cout << "start iteration... " << endl;
		mod.solve_selfcon(sim.mode_calc);
		mod.find_transitions(simP);
	//
		// cout << "solving approx... " << endl;
		if ((sim.mode_stats == 2) || (sim.mode_stats == 3))
		{
			// mod_approx = mod;
			mod_approx.solve_selfcon(1);
			mod_approx.find_transitions(simP);

			compare_approx(modP,modP_approx);
		}
		// cout << "done" << endl;

		if (sim.mode_stats == 4 && !modP->simulation.in_inc)
		{
			mod.integrate_information();
		}

		// cout << "storing ...";
		sim.store_results(modP,modP_approx);
		// cout << "done" << endl;
	};
	cout << "computation done!" << endl;
	//
	// // if (sim.steps == 1)
	// 	// write_theory(fileOut);
	// // else {
    write_results(fileOut,simP,modP,modP_approx);
    // // }
	return 0;
}
