
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
#include "./computations/comFunctions.cpp"
#include "./io/readwrite.cpp"

int main(int argc, char** argv)
{

	static plog::ConsoleAppender<plog::TxtFormatter> consoleAppender;
	static plog::RollingFileAppender<plog::CsvFormatter> fileAppender("MultiAppender.csv", 8000, 3);
	static plog::Severity sev = plog::warning; 

	// spdlog::set_level(spdlog::level::info);
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

	if ((argc < 4) or (argc > 5))
	{
		PLOG_ERROR << "Please specify (only) the input and output file names!" << endl;
		return 1;
	}


	string fileModel = argv[1];
    string fileCom = argv[2];
	string fileOut = argv[3];

    Model mod, *modP = &mod;
    Computation com, *comP = &com;
    Measures mes, *mesP = &mes;

    read_model(fileModel,modP);
    read_computation(fileCom,comP);

    mod.set_weights();
    mod.solve_selfcon(0);
    // mod.find_transitions(simP);
	mod.get_max_prob();

	mes.rates.resize(mod.layer[0].nPop);
	mes.rates_T.resize(mod.layer[0].nPop);

	// initiate random number generator
	gsl_rng_env_setup();
    const gsl_rng_type *TYPE = gsl_rng_mt19937;
    gsl_rng *rng = gsl_rng_alloc(TYPE);
	gsl_rng_set(rng, 1 + com.seed);

	// cout << "drawing firing rates from theoretical distribution..." << endl;
	for (unsigned p=0; p<mod.layer[0].nPop; p++) {
		PLOG_DEBUG << "population " << p << ": drawing firing rates from theoretical distribution..." << endl;
		Population_Simulation *popSimP = &mod.layer[0].population[p].simulation;
		mes.rates[p].resize(com.N);
		mes.rates_T[p].resize(com.N);

		for (unsigned n=0; n<com.N; n++) {

			PLOG_DEBUG << "population " << p << ", neuron " << n << ": drawing firing rates from theoretical distribution..." << endl;
			mes.rates[p][n].resize(com.draw_from_theory);
			mes.rates_T[p][n].resize(com.draw_from_theory);

			for (unsigned k=0; k<com.draw_from_theory; k++) {
				mes.rates[p][n][k] = popSimP->draw_rate(rng);
				mes.rates_T[p][n][k].resize(com.draw_finite_time);

				for (unsigned t=0; t<com.draw_finite_time; t++) {
					mes.rates_T[p][n][k][t] = popSimP->draw_sample(rng,mes.rates[p][n][k],com.T);
					PLOG_DEBUG << "neuron n with rate " << mes.rates[p][n][k] << ", drew rate: " << mes.rates_T[p][n][k][t] << endl;
				}
			}
		}
	}

    // res.p_bayes_est.resize(com.draw_from_theory);

	//
	//
	// 	cout << "\r\tk=" << com.k+1 << "/" << com.draw_from_theory << flush;
	//
    //     // obtain density estimate from that
    //     // res.p_est_inf = get_density_estimate(res.N_AP, com, "poisson");
	//
    //     // now, obtain rate of neuron n, measured in finite time intervall T
    //     // if (com.draw_finite_time > 0)
    //     // {
    //     // resP->p_bayes_est[comP->k].resize(comP->draw_finite_time);
	//
    //     // cout << "drawing samples of finite time measurements with T=" << com.T << ", and calculating bayes estimate of firing rates using the '" << com.prior << "' prior. "<< endl;
    //     for (com.j=0; com.j<com.draw_finite_time; ++com.j) {
    //         com.draw_samples(mesP);
	//
    //         // obtain density estimate from that
    //         // bayesian_estimate_theory(modP, comP, resP);
    //     }
    // }

	write_measures(fileOut, comP, modP, mesP);

    return 0;
}
