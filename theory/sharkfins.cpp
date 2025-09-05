/*
	master-code to run all setup, iteration and solving functions from

	created by Alexander Schmidt
	last changed on November 1st 2023

	TODO:
		* set up proper logging to include debugging outputs when called accordingly
*/ 


#include <plog/Log.h>
#include <plog/Init.h>

#include <plog/Formatters/CsvFormatter.h> 
#include <plog/Appenders/RollingFileAppender.h> 
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ConsoleAppender.h> 

#include "./computations/structures.h"
#include "./computations/functions.h"
#include "./computations/functions.cpp"
#include "./computations/modFunctions.cpp"
#include "./computations/simFunctions.cpp"
#include "./io/readwrite.cpp"

int main(int argc, char** argv)
{
	// check, if proper number of arguments was provided

	static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
	static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("MultiAppender.csv", 8000, 3);
	static plog::Severity sev = plog::warning; 

	if (argc==5) {
		// set up logger to write to console
		switch (atoi(argv[4])) {
			case 5:	//'debug':
				sev = plog::debug;
				break;
			case 4:	//'info':
				sev = plog::info;
				break;
			case 3: //'warning':
				sev = plog::warning;
				break;
			case 2:	//'error':
				sev = plog::error;
				break;
		}	
	}
	
	plog::init(sev, &consoleAppender);//.addAppender(&fileAppender);

	if ((argc < 3) or (argc > 5))
	// if (argc != 4)
	{
		PLOG_ERROR << "Please specify (only) the input and output file names!" << endl;
		return 1;
	}

	// obtain paths of input files
	unsigned file_idx = 1;
	string fileModel = argv[file_idx++];
	string fileSim;
	if (argc>3) {
		fileSim = argv[file_idx++];
	}
	string fileOut = argv[file_idx++];
	
	
	// prepare and read paramaters from input files
	PLOG_DEBUG << "Initiating model arrays " << endl;
    Model mod, *modP = &mod;
    read_model(fileModel,modP);

	Model mod_approx, *modP_approx = &mod_approx;
    read_model(fileModel,modP_approx);
	PLOG_DEBUG << "Initiating model arrays ... done" << endl;
	
	Simulation sim, *simP = &sim;
	if (argc>3) {
		read_simulation(fileSim,simP);
		sim.initialize(modP);
		sim.initialize(modP_approx);
	};

	// initialize arrays and functions
	PLOG_DEBUG << "Initiating results arrays " << endl;
    initiate_results(modP,simP);
    initiate_results(modP_approx,simP);

	PLOG_DEBUG << "Initiating results arrays ... done" << endl;

	// cout << "mode selfcon: " << sim.mode_selfcon << endl;

	// run iteration to solve model for each parameter set
	// PLOG_DEBUG << "calculate solutions for mode stats = " << sim.mode_stats;
	
	do {
		// PLOG_DEBUG << "start iteration... ";

		mod.set_parameters();
		mod_approx.set_parameters();

		// PLOG_DEBUG << "solving model... ";
		
		if (sim.mode_selfcon==0)
		mod.solve_selfcon(sim.mode_calc);
		else 
		mod.solve_selfcon_from_currents(sim.mode_calc);
		
		// PLOG_DEBUG << "finding transition... ";
		mod.find_transitions(simP);

		if ((sim.mode_stats == 2) || (sim.mode_stats == 3)) {
			// PLOG_DEBUG << "solving approximate model... ";
			mod_approx.solve_selfcon(1);
			// PLOG_DEBUG << "finding transitions in approx model";
			mod_approx.find_transitions(simP);

			// PLOG_DEBUG << "comparing approx model";
			compare_approx(modP,modP_approx);
			// PLOG_DEBUG << "done";
		}

		if (sim.mode_stats == 4 && !modP->simulation.in_inc)
			mod.integrate_information();

		PLOG_DEBUG << "storing ...";
		sim.store_results(modP,modP_approx);
		PLOG_DEBUG << "done";
	} while (sim.run_iteration(modP,modP_approx));

	PLOG_DEBUG << "computation done!";

    write_results(fileOut,simP,modP,modP_approx);
	return 0;
}
