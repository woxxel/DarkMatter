// #include <vector>
// #include <math.h>
// #include <stdarg.h>
// #include <time.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_multiroots.h>
#include <random>
#include <gsl/gsl_randist.h>
#include "gsl/gsl_rng.h"
// #include <gsl/gsl_integration.h>

#include "functions.h"
// #include "structures.h"

double poisson_distr(int k, double lambda)
{
        return exp(k*log(lambda)-lambda-lgamma(k+1.0));
}

double Population_Simulation::draw_rate(gsl_rng *rng)
{
    // initiate variables
    double rate;
    int i=0;    // track iterations to break off in case of inability to find something

    // rejection sampling algorithm
    while (true) {
        rate = gsl_rng_uniform_pos(rng)*rate_max; // generate random rate
        // and accept sample according to theoretically predicted probability
        if (gsl_ran_flat(rng, 0,max_prob) < distribution_exact(rate))
            return rate;
            // mesP->rates[k][n] = rate_tmp;
            // resP->N_AP[comP->k][n] = int(rate_inf_tmp*comP->T);
        if (i++ > 100000) return NAN;
    }
}

double Population_Simulation::draw_sample(gsl_rng *rng, double rate, double T)
{
    // int N_AP = gsl_ran_poisson(rng,rate*T);
    return gsl_ran_poisson(rng,rate*T)/T;     // generating a random number of spikes from this
}

// vector<double> get_density_estimate(vector<int> data, computation sim, string kernel)
// {
//         vector<double> p_est(simP->AP_max,0);
//         cout << "getting density estimation..." << endl;
//     //     vector<vector<double> > p_est_single (simP->N);
//         for (int n=0; n<simP->N; ++n)
//         {
//                 if (kernel.compare("poisson") == 0)
//                         for (int i=0; i<simP->AP_max; ++i)
//                                 p_est[i] += poisson_distr(data[n],i)/(simP->N/simP->T);
//                 else
//                         cout << "kernel '" << kernel << "' not yet implemented" << endl;
//         }
//     //     for (int i=0; i<N_max; ++i)
//     //         cout << "prob(i=" << i << "): " << p_est[i] << endl;
//         return p_est;
// }


// void bayesian_estimate(model *modP, computation *comP, results *resP)
// {
//         // wanna calculate: P(lambda|k) = P(k|lambda)*p(lambda)/p(k)
//         // if no prior beliefs exist, p(lambda) = p(k) -> P(lambda|k) = P(k|lambda)
//         // however, if we take p(lambda) to be firing rate distribution from theory, we obtain the following:
//
// //     P_k_lambda  = poisson_distr(k,nu*T);
// //     p_lambda    = modP->distribution_exact(nu,0);
// //     p_k = int(P_k_lambda*p_lambda) d lambda
//
//         resP->p_bayes_est[comP->k][comP->j].resize(simP->steps,0);
//         for (int n=0; n<comP->N; ++n)
//         {
// //                 cout << "APs for neuron " << n << ": " << k << endl;
//                 double p_k = 0;
//
//                 for (int i=1; i < simP->steps; ++i)    // "integration" over lambda    //should also work with integration (faster?)
//                         p_k += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior);
//
//                 resP->p_bayes_est[comP->k][comP->j][0] = 0;
//                 for (int i=1; i < simP->steps; ++i)    // construction of posterior
//                         resP->p_bayes_est[comP->k][comP->j][i] += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior)/p_k / (resP->d_nu*comP->N);
//         }
// }
// //
// //
// double bayes_est_prior(double nu, model *modP, string prior)
// {
//         if (not prior.compare("mean_field"))
//             return modP->distribution_exact(nu,0);
//         else if (not prior.compare("plain"))
//             return 1;
//         else
//         {
//             cout << "Specify a prior!" << endl;
//             exit(0);
//         }
// }


vector<double> get_cdf(vector<double> p, int steps, double d_nu)     // compute cdf
{
        vector<double> cdf(steps,0);
        for(int i=1; i<=steps; ++i)
                cdf[i] = cdf[i-1] + p[i]*d_nu;

        return cdf;
}

// void post_process(computation *comP, model *modP, model *mod_approxP, simulation *simP, results *resP)
// {
// //         resP->KS.resize(comP->draw_from_theory);
// //         resP->KL_entropy_single.resize(simP->steps);
//
// //         cout << "gamma: " << modP->paras.gamma[0] << "," << mod_approxP->paras.gamma[0] << endl;
// //         cout << "delta: " << modP->paras.delta[0] << "," << mod_approxP->paras.delta[0] << endl;
//
// //         double  KL_tmp = KL_divergence(0,modP->paras.rate_max[0],modP->paras,mod_approxP->paras);
// //         cout << "KL = " << KL_tmp << endl;
// //         for (int i=1; i<simP->steps; ++i)
// //         {
// //                 resP->KL_entropy_single[i] = KL_divergence((i-1)*resP->d_nu,i*resP->d_nu,modP->paras,mod_approxP->paras);
// //                 cout << "KL[" << i << "] = " << resP->KL_entropy_single[i] << endl;
// //         }
//
//         for (int k=0; k<comP->draw_from_theory; ++k)
//         {
//                 resP->KS[k].resize(comP->draw_finite_time,0);
// //                 resP->KL_entropy[k].resize(comP->draw_finite_time,0);
//
//                 for (int j=0; j<comP->draw_finite_time; ++j)
//                 {
//                         vector<double> cdf_tmp = get_cdf(resP->p_bayes_est[k][j],simP->steps,resP->d_nu);
//
//                         double KS_tmp = 0;
//                         for (unsigned i=1; i<simP->steps; ++i)
//                         {
//                                 KS_tmp = abs(cdf_tmp[i] - resP->cdf_theory[i]);
// //                                 cout << "KS tmp: " << KS_tmp << ",\tKS old: " << resP->KS[k][j] << endl;
//                                 resP->KS[k][j] = (KS_tmp > resP->KS[k][j]) ? KS_tmp : resP->KS[k][j];
//
// //                                 cout << "p[" << i << "]: " << resP->p_exact[i] << " vs " << resP->p_bayes_est[k][j][i] << ", log2: " << log2(resP->p_exact[i]/resP->p_bayes_est[k][j][i]) << endl;
// //                                 if(resP->p_bayes_est[k][j][i]>pow(10,-3) && resP->p_exact[i]>pow(10,-3))
// //                                     resP->KL_entropy[k][j] += resP->p_exact[i] * resP->d_nu * log2(resP->p_exact[i]/resP->p_bayes_est[k][j][i]);
//                         }
// //                         cout << "KS test: " << resP->KS[k][j] << endl;
//                 }
//         }
//
// }
