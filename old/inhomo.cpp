#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdio.h>
#include <math.h>

#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "functions.cpp"
#include "readwrite.cpp"

// using namespace std;

int main(int argc, char** argv)
{	
	if ((argc < 4) or (argc > 4))
	{
		cout << "Please specify (only) the input and output file name!" << endl;
		return 1;
	}
	
	string fileModel = argv[1];
        string fileSim = argv[2];
	string fileOut = argv[3];
        
        simulation sim, *simP = &sim;
        model mod, *modP = &mod;
        results res;//, *resP = &res;
        
        read_model(fileModel,modP);
        read_simulation(fileSim,simP);
        
        //! initialize paras as vector, containing the set of network parameters for the single iterations, each
        //! therefore, get size first!
        
        double val = -1;
        
        // a little dirty...
        double border_inc[sim.steps] = {NAN};//{[0 ... sim.steps] = NAN};
        for (unsigned i=0;i<sim.steps;++i)
            border_inc[i] = NAN;

	double gamma[sim.nSz][sim.alpha_0Sz][sim.tau_GSz][sim.rateWntSz][mod.paras.population];
	double *gammaP = &gamma[0][0][0][0][0];
	double chi[sim.nSz][sim.alpha_0Sz][sim.tau_GSz][sim.rateWntSz][mod.paras.population];
	double *chiP = &chi[0][0][0][0][0];
	
        for (sim.alpha_0_iter = 0; sim.alpha_0_iter < sim.alpha_0Sz; sim.alpha_0_iter++)
        {
                mod.paras.alpha_0 = sim.alpha_0[sim.alpha_0_iter];
                
                for (sim.rateWnt_iter = 0; sim.rateWnt_iter < sim.rateWntSz; sim.rateWnt_iter++)
                {
                        mod.paras.rateWnt = sim.rateWnt[sim.rateWnt_iter];
                        
// 			vector<bool> found_inconsistent (mod.paras.population,false);
                        vector<bool> found_no_peak (mod.paras.population,false);
                        
// 			cout << "nu: " << sim.rateWnt[sim.rateWnt_iter] << " ,\t alpha_0: " << sim.alpha_0[sim.alpha_0_iter] << endl;
                        
                        mod.solve_selfcon("exact");
                        
                        //! post processing
                        
                        // search all populations for inconsistent states. If found, fill following gamma with "-1"
                        for (unsigned p = 0; p < mod.paras.population; p++)
                        {
                                if (mod.q_border1(p) || mod.q_border2(p))// && (found_inconsistent[p] == 0))
                                {
//                                                         cout << "found inconsistent" << endl;
                                        cout << "need to get first value of inconsistence here" << endl;
                                        // check, which one is the x-axis
                                        /*double val = -1;*///nan("");
// 							if (sim.rateWntSz > 1)
//                                                         {
//                                                                 border_inc[sim.rateWnt_iter] = sim.y_val();
// 								for (unsigned f = sim.rateWnt_iter; f < sim.rateWntSz; f++)
// 									gamma[sim.n_iter][sim.alpha_0_iter][sim.tau_G_iter][f][p] = val;
//                                                         }
// 							else if (sim.tau_GSz > 1)
//                                                         {
//                                                                 border_inc[sim.tau_G_iter] = sim.y_val();
// 								for (unsigned f = sim.tau_G_iter; f < sim.tau_GSz; f++)
// 									gamma[sim.n_iter][sim.alpha_0_iter][f][sim.rateWnt_iter][p] = val;
//                                                         }
// 							else if (sim.alpha_0Sz > 1)
//                                                         {
//                                                                 border_inc[sim.alpha_0_iter] = sim.y_val();
// 								for (unsigned f = sim.alpha_0_iter; f < sim.alpha_0Sz; f++)
// 									gamma[sim.n_iter][f][sim.tau_G_iter][sim.rateWnt_iter][p] = val;
//                                                         }
//                                                         else
//                                                         {
//                                                                 cout << "no x-axis could be found - quitting..." << endl;
//                                                                 exit(0);
//                                                         }
                                }
                                else
                                        gamma[sim.n_iter][sim.alpha_0_iter][sim.tau_G_iter][sim.rateWnt_iter][p] = mod.paras.gamma[p];
                        }
// 					if (mod.inconsistent())
//                                         {
// //                                                 cout << "found network to be inconsistent" << endl;
// 						break;
// 					}
                        
                        for (unsigned p = 0; p < mod.paras.population; p++)
                        {
                                if(mod.no_peak(p))
                                {
                                        cout << "need to get first value of no peak here" << endl;
                                }
                                        found_no_peak[p] = 1;
// 					cout << "no peak" << endl;
// 					for (unsigned k = j; k < sim.rateWntSz; k++)
                                        gamma[sim.n_iter][sim.alpha_0_iter][sim.tau_G_iter][sim.rateWnt_iter][p] = -2;			// entries -1 mean: inconsistent
                                }
// 						else
// 						{
// 							double nu_peak_log_I = mod.nu_peak_log_full(mod.paras.gamma[p],mod.paras.delta[p],mod.paras.rate_max[p]);
                                        chi[sim.n_iter][sim.alpha_0_iter][sim.tau_G_iter][sim.rateWnt_iter][p] = mod.paras.chi[p];
// 						}
                        }
                }
// 			}
        }
// 	}
	
// 	for (unsigned i=0;i<sim.steps;++i)
//             cout << "border: " << border_inc[i] << endl;
	
// 	gsl_multiroot_fsolver_free (s);
// 	gsl_vector_free (q);
	write_sharks(fileOut,sim,mod,gammaP,chiP);
	return 0;
}

// void print_state (size_t iter, gsl_multiroot_fsolver * s)
// {
//   printf ("iter = %3u x = % .3f % .3f " "f(x) = % .3e % .3e\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
// }