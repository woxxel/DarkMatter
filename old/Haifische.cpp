#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdio.h>

#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "functions.cpp"
#include "read_haifische.cpp"

// using namespace std;

int main(int argc, char** argv)
{
	if ((argc < 3) or (argc > 3))
	{
		cout << "Please specify (only) the input and output file name!" << endl;
		return 1;
	}

	string fileModel = argv[1];
	string fileOut = argv[2];

        simulation sim, *simP = &sim;
        model mod, *modP = &mod;
        results res, *resP = &res;

	int status;

        read_model(fileModel,modP);
        read_sim(fileSim,simP);

        //! initialize paras as vector, containing the set of network parameters for the single iterations, each
        //! therefore, get size first!

	double I_balance, delta;

	double gamma[sim.n.size()][sim.alpha_0.size()][sim.tau_I.size()][sim.rateWnt.size()][2];
	double *gammaP = &gamma[0][0][0][0][0];
	double chi[sim.n.size()][sim.alpha_0.size()][sim.tau_I.size()][sim.rateWnt.size()][2];
	double *chiP = &chi[0][0][0][0][0];

	for (unsigned n_iter = 0; n_iter < sim.n.size(); n_iter++)
	{
		mod.paras.n = sim.n[n_iter];

		for (unsigned alpha_iter = 0; alpha_iter < sim.alpha_0.size(); alpha_iter++)
		{
			for (unsigned p = 0; p < 2; p++)
				net.paras.alpha_0[p] = sim.alpha_0[alpha_iter];		// set rates
// 			cout << "alpha: " << net.paras.alpha_0[0] << endl;
			for (unsigned tau_I_iter = 0; tau_I_iter < sim.tau_I.size(); tau_I_iter++)
			{
				net.paras.tau_I = sim.tau_I[tau_I_iter];
// 				cout << "tau: " << net.paras.tau_I << endl;

				for (unsigned rate_iter = 0; rate_iter < sim.rateWnt.size(); rate_iter++)
				{
					vector<int> found_inconsistent (2,0);
					vector<int> found_no_peak (2,0);

// 					gsl_vector *q_guess = gsl_vector_alloc (2);
// 					for (unsigned p = 0; p < 2; p++)
// 					{
// 						net.paras.rateWnt[p] = sim.rateWnt[rate_iter];				// set threshold heterogeneity
// 						gsl_vector_set(q_guess, p, gsl_pow_2(net.paras.rateWnt[p]));	// first guess
// 					}
//
// // 					cout << "nu: " << sim.rateWnt[rate_iter] << ", alpha_0: " << sim.alpha_0[alpha_iter] << endl;
//
// 					gsl_multiroot_function F;
// 					F.f = &selfconsistency_f;
// 					F.n = 2;
// 					F.params = &net.paras;
//
//
// 					net.get_sigma_V();
// 					gsl_multiroot_fsolver_set (s, &F, q_guess);
//
// 					size_t iter = 0;
//
// 					do
// 					{
// 						iter++;
// 						status = gsl_multiroot_fsolver_iterate (s);
//
// 						if (status) break;	// check if solver is stuck
//
// 						status = gsl_multiroot_test_residual (s->f, 1e-7);
// 					}
// 					while (status == GSL_CONTINUE && iter < 100);
		// 			print_state (iter, s);
		// 			printf ("status = %s\n", gsl_strerror (status)); // implement here breaking when search fails

// 					net.paras.q[0] = gsl_vector_get(s->x,0);
// 					net.paras.q[1] = gsl_vector_get(s->x,1);

// 					cout << "q: " << net.paras.q[0] << endl;
					//! post processing

		// 			cout << "solved: q_I=" << net.paras.q[0] << ", q_E=" << net.paras.q[1] << endl;

					for (unsigned p = 0; p < 2; p++)
					{
						net.paras.alpha[p] = sqrt( gsl_pow_2(net.paras.J_I[p]) * net.paras.q[0] + gsl_pow_2(net.paras.J_E[p]) * net.paras.kappa * net.paras.q[1] + gsl_pow_2(net.paras.alpha_0[p]) );

						double q_border1 = net.paras.rateWnt[p]/net.paras.rate_max[p] * sqrt(gsl_pow_2(net.paras.alpha[p]) + gsl_pow_2(net.paras.sigma_V[p]))/net.paras.sigma_V[p];
						double q_border2 = net.paras.q[p]/gsl_pow_2(net.paras.rate_max[p]) * sqrt(2 * gsl_pow_2(net.paras.alpha[p]) + gsl_pow_2(net.paras.sigma_V[p]))/net.paras.sigma_V[p];

						if (((1 < q_border1) || (1 < q_border2)) && (found_inconsistent[p] == 0))
						{
							found_inconsistent[p] = 1;

							double val = -1;//nan("");
							if (sim.rateWnt.size() > 1)
								for (unsigned f = rate_iter; f < sim.rateWnt.size(); f++)
									gamma[n_iter][alpha_iter][tau_I_iter][f][p] = val;
							else if (sim.tau_I.size() > 1)
								for (unsigned f = tau_I_iter; f < sim.tau_I.size(); f++)
									gamma[n_iter][alpha_iter][f][rate_iter][p] = val;
							else if (sim.alpha_0.size() > 1)
								for (unsigned f = alpha_iter; f < sim.alpha_0.size(); f++)
									gamma[n_iter][f][tau_I_iter][rate_iter][p] = val;
						}
					}
					if ((gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][0] == -1) and (gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][1] == -1))
						break;

					for (unsigned p = 0; p < 2; p++)
					{
						I_balance = sqrt(I_squared_nu(net.paras.alpha[p],net.paras.sigma_V[p],net.paras.rateWnt[p],net.paras.rate_max[p]));
						delta = I_balance/net.paras.alpha[p];

						gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][p] = net.paras.sigma_V[p]/net.paras.alpha[p];

						if(gsl_pow_2(gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][p])*gsl_pow_2(delta)-4*(gsl_pow_2(gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][p]) - 1) < 0)
						{
							found_no_peak[p] = 1;
		// 					cout << "no peak" << endl;
		// 					for (unsigned k = j; k < sim.rateWnt.size(); k++)
							gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][p] = -2;			// entries -1 mean: inconsistent
						}
						else
						{
							double nu_peak_log_I = nu_peak_log_full(gamma[n_iter][alpha_iter][tau_I_iter][rate_iter][p],delta,net.paras.rate_max[p]);
							chi[n_iter][alpha_iter][tau_I_iter][rate_iter][p] = -log10(exp(1)) * nu_peak_log_I + log10(net.paras.rateWnt[0]);
						}
					}
				}
			}
		}
	}
	gsl_multiroot_fsolver_free (s);
// 	gsl_vector_free (q);
	write_results(fileOut,sim,gammaP,chiP);
	return 0;
}

// void print_state (size_t iter, gsl_multiroot_fsolver * s)
// {
//   printf ("iter = %3u x = % .3f % .3f " "f(x) = % .3e % .3e\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
// }
