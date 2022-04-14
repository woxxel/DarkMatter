#include "./computations/structures.h"
#include "./computations/functions.h"
#include "./computations/functions.cpp"
#include "./computations/modFunctions.cpp"
#include "./computations/simFunctions.cpp"
#include "./computations/comFunctions.cpp"
#include "./io/readwrite.cpp"

int main(int argc, char** argv)
{

	// spdlog::set_level(spdlog::level::info);
	if ((argc < 4) or (argc > 4))
	{
		cout << "Please specify (only) the input and output file names!" << endl;
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

	mes.rates.resize(com.N);
	mes.rates_T.resize(com.N);

	// initiate random number generator
	gsl_rng_env_setup();
    const gsl_rng_type *TYPE = gsl_rng_mt19937;
    gsl_rng *rng = gsl_rng_alloc(TYPE);
	gsl_rng_set(rng, 1 + com.seed);

	// cout << "drawing firing rates from theoretical distribution..." << endl;
	for (unsigned n=0; n<com.N; n++) {
		mes.rates[n].resize(com.draw_from_theory);
		mes.rates_T[n].resize(com.draw_from_theory);

		for (unsigned k=0; k<com.draw_from_theory; k++) {
			mes.rates[n][k] = com.draw_rate(rng,modP);
			mes.rates_T[n][k].resize(com.draw_finite_time);

			for (unsigned t=0; t<com.draw_finite_time; t++) {
				mes.rates_T[n][k][t] = com.draw_sample(rng,mes.rates[n][k],com.T);
				// cout << "neuron n with rate " << mes.rates[n][k] << ", drew rate: " << mes.rates_T[n][k][t] << endl;
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

	write_measures(fileOut,comP,mesP);

    return 0;
}
