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
    string fileSim = argv[2];
	string fileOut = argv[3];

    Model mod, *modP = &mod;
    Computation com, *comP = &com;
    Measures mes, *mesP = &mes;

    read_model(fileModel,modP);
    read_computation(fileSim,comP);

    mod.set_weights();
    mod.solve_selfcon(0);
    // mod.find_transitions(simP);

    if (com.draw_from_theory > 0)
    {
        mes.rates.resize(com.draw_from_theory);
        mes.rates_T.resize(com.draw_from_theory);
        mes.N_AP.resize(com.draw_from_theory);
        // res.p_bayes_est.resize(com.draw_from_theory);

        cout << "drawing firing rates from theoretical distribution..." << endl;
        for (com.k=0; com.k<com.draw_from_theory; ++com.k)
        {
            cout << "\r\tk=" << com.k+1 << "/" << com.draw_from_theory << flush;
            com.draw_rates(modP,mesP);

            // obtain density estimate from that
            // res.p_est_inf = get_density_estimate(res.N_AP, com, "poisson");

            // now, obtain rate of neuron n, measured in finite time intervall T
            if (com.draw_finite_time > 0)
            {
                mes.rates_T[com.k].resize(com.draw_finite_time);
                mes.N_AP[com.k].resize(com.draw_finite_time);
                // resP->p_bayes_est[comP->k].resize(comP->draw_finite_time);

                // cout << "drawing samples of finite time measurements with T=" << com.T << ", and calculating bayes estimate of firing rates using the '" << com.prior << "' prior. "<< endl;
                for (com.j=0; com.j<com.draw_finite_time; ++com.j)
                {
                            // cout << "\tj=" << com.j+1 <<"/" << com.draw_finite_time << endl;
                    com.draw_samples(mesP);

                    // obtain density estimate from that
                    // bayesian_estimate_theory(modP, comP, resP);
                }
            }
        }
        cout << endl;
    }

    return 0;
}
