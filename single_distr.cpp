#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdio.h>
#include <random>
#include <math.h>

#include "functions.cpp"
#include "bayes_estimate.cpp"
#include "readwrite.cpp"



int main(int argc, char** argv)
{
        measures mes, *mesP = &mes;
        model mod, *modP = &mod;
        model mod_approx, *mod_approxP = &mod_approx;
        simulation sim, *simP = &sim;
        computation com, *comP = &com;
        results res, *resP = &res;
        
        string fileMeasures, fileModel, fileSim, fileComp, fileOut;

        if (atoi(argv[1]) == 0)
        {
            cout << "Reading model and computation parameters to generate measurements..." << endl;
            fileModel = argv[2];
            fileSim = argv[3];
            fileComp = argv[4];
            fileOut = argv[5];
            
            read_model(fileModel,modP);
            read_model(fileModel,mod_approxP);
            read_simulation(fileSim,simP);
            read_computation(fileComp,comP);
        }
        else if (atoi(argv[1]) == 1)
        {
            cout << "Reading data from measurements and computation parameters..." << endl;
            fileMeasures = argv[2];
            fileComp = argv[3];
            fileOut = argv[4];
            
            read_measures(fileMeasures,mesP);
            read_computation(fileComp,comP);
        }
        else
        {
            cout << "Please specify (only) the input and output file names!" << endl;
            return 1;
        }
        
        
        if (com.p_theory > 0)
        {
            mod.set_weights();
            // get proper d_nu to not waste too much ressources
            if (com.draw_from_theory > 0)
            {
                    res.steps = 1000;
                    res.factor = 0.5;
                    while (res.steps >= 1000)
                    {
                        res.factor *= 2;
                        res.steps = sim.rateWnt[0]*com.border*com.T/res.factor;
        //                 cout << "steps: " << res.steps << endl;
                    }
            }
            else
                    res.steps = com.T;
            for (int p = 0; p < mod.paras.Npop; p++)
            {
                    mod.paras.rate[p] = sim.rateWnt[0];
                    mod.paras.tau_G = sim.tau_G[0];
                    mod.paras.alpha_0[p] = sim.alpha_0[0];
                    
            }
            res.d_nu = sim.rateWnt[0]*com.border/res.steps;
            
            cout << "steps: " << res.steps << ", d_nu=" << res.d_nu << endl;
            
            //! obtain analytical firing rate distribution (Hiemeyer)
            cout << "J: " << mod.paras.J_I[0] << ",\t rate: " << sim.rateWnt[0] << ",\t tau_G: " << sim.tau_G[0] << ",\t tau_M: " << mod.paras.tau_M << endl;
            
            mod.solve_selfcon(0);
            cout << "exact:  \t" << "alpha: " << mod.paras.alpha[0] << " ,\t sigma: " << mod.paras.sigma_V[0] << " ,\t gamma': " << mod.paras.gamma[0] << endl;
            
//             mod_approx.solve_selfcon(1);
//             cout << "approx: \t" << "alpha': " << mod_approx.paras.alpha[0] << " ,\t sigma': " << mod_approx.paras.sigma_V[0] << ",\t gamma': " << mod_approx.paras.gamma[0] << endl;
            
            res.alpha_single = mod.paras.alpha[0];
            res.sigma_V_single = mod.paras.sigma_V[0];
            res.I_single = -sqrt(I_squared_nu(res.alpha_single,res.sigma_V_single,sim.rateWnt[0],mod.paras.rate_max[0]));
            res.chi_single = mod.paras.chi[0];
//             mod.solve_selfcon("approx");
            
//             // write parameter to results-struct
//             res.delta.resize(mod.paras.population);
//             res.gamma.resize(mod.paras.population);
//             res.chi.resize(mod.paras.population);
            
//             for (unsigned p = 0; p < mod.paras.population; p++)
//             {
//                     res.delta[p] = mod.paras.delta[p];
//                     res.gamma[p] = mod.paras.gamma[p];
//                     res.chi[p] = mod.paras.chi[p];
//             }
            
            double cosh_tmp = cosh(mod.paras.gamma[0]*mod.paras.delta[0]*sqrt(-2*log(res.d_nu/mod.paras.rate_max[0]))); // check if pdf can be calculated
            if (!isinf(cosh_tmp) and (com.p_theory==2))
            {
                    // write function values and x-axis values into results array
                    cout << "Calculate pdf and cdf from selfconsistent solution...";
                    
                    res.p_range.resize(res.steps);
                    res.p_exact.resize(res.steps);

                    res.p_exact[0] = 0;
                    for (int i=1; i<res.steps; ++i)
                    {
                            res.p_range[i] = i*res.d_nu;    // should be written in some seperate function
                            res.p_exact[i] = mod.distribution_exact(res.p_range[i],0);
                            res.max_prob = max(res.p_exact[i],res.max_prob);      // obtain maximum value of probability density (with bins of width=rate_max/N_max)
                            
//                             cout << "i: " << i << " ,\t nu=" << res.p_range[i] << " ,\t p(nu)=" << res.p_exact[i] << " ,\t max=" << res.max_prob << endl;
                    }
                    
//                     res.p_approx.resize(res.steps);
//                     res.p_approx[0] = 0;
//                     for (int i=1; i<res.steps; ++i)
//                             res.p_approx[i] = mod_approx.distribution_exact(res.p_range[i],0);
                    
                    cout << "done!" << endl;
                                        
                    // transform pdf into a histogram by integrating over d_nu (put to function)
                    if (com.p_theory_hist > 0)
                    {
                            cout << "Transform pdf into a histogram with expected rates drawn from a population of N neurons...";
                            res.p_hist.resize(com.n_bin);
                            int k=0;
                            int factor = 1/res.d_nu;
                            
                            for (int i=0; i<res.steps; ++i)
                                    if (i%factor == 0)
                                    {
                                            res.p_hist[k] = round(pdf2hist(res.d_nu*k,res.d_nu*(k+1),mod.paras)*com.N); // is fully defined, as gamma and delta are written in paras
                                            ++k;
                                    }
                            cout << "done!" << endl;
                    }
                    
                    cout << "Calculating cdf...";
                    res.cdf_theory.resize(res.steps,0);
                    for(int i=1; i<=res.steps; ++i)
                    {
                        res.cdf_theory[i] = pdf2hist((i-1)*res.d_nu,i*res.d_nu,mod.paras)+res.cdf_theory[i-1];
                        cout << "cdf[" << i << "] = " << res.cdf_theory[i] << endl;
                    }
                    cout << "done!" << endl;
            }
//             else
//                     cout << "pdf cannot be calculated, as cosh exceeds numerical range" << endl;
            
            //! obtain rates of N neurons sampled from theoretically predicted trajectory (by rejection sampling)
            if (com.draw_from_theory > 0)
            {
                res.rate_inf.resize(com.draw_from_theory);
                res.p_bayes_est.resize(com.draw_from_theory);

                if (com.draw_finite_time > 0)
                {
                    res.rate_T.resize(com.draw_from_theory);
                    res.N_AP.resize(com.draw_from_theory);
                }
                
                cout << "drawing firing rates from theoretical distribution..." << endl;
                for (com.k=0; com.k<com.draw_from_theory; ++com.k)
                {
//                     res.rate_inf.push_back(vector<double>(com.N,0));
//                     cout << "\r\tk=" << com.k+1 << "/" << com.draw_from_theory << flush;
                    draw_rates(modP,comP,resP);
                    
                    // obtain density estimate from that
//                     res.p_est_inf = get_density_estimate(res.N_AP, com, "poisson");
                    
                    // now, obtain rate of neuron n, measured in finite time intervall T
                    if (com.draw_finite_time > 0)
                    {
                        res.rate_T[com.k].resize(com.draw_finite_time);
                        res.N_AP[com.k].resize(com.draw_finite_time);
                        resP->p_bayes_est[comP->k].resize(comP->draw_finite_time);
                        
//                         cout << "drawing samples of finite time measurements with T=" << com.T << ", and calculating bayes estimate of firing rates using the '" << com.prior << "' prior. "<< endl;
                        for (com.j=0; com.j<com.draw_finite_time; ++com.j)
                        {
//                             cout << "\tj=" << com.j+1 <<"/" << com.draw_finite_time << endl;
                            draw_samples(comP,resP);
                            
                            // obtain density estimate from that
                            bayesian_estimate_theory(modP, comP, resP);
                            
//                             post_process_measure(resP);         // calculate KL-entropy and KS-test and bootstrap?!
//                             res.p_est_T = get_density_estimate(res.N_AP, com, "poisson");
                        }
                    }
                }
                cout << endl;
            }
//             cout << "post process" << endl;
            
            post_process(comP,modP,mod_approxP,resP);
//             cout << "writing stuff" << endl;
            write_theory(fileOut, resP);
        }
        else if(com.process_data > 0)
        {
            int AP_max = 0;
            
            for (int n=0; n<mes.N; ++n)
                AP_max = (mes.N_AP[n] > AP_max) ? mes.N_AP[n] : AP_max;
            
            res.steps = 2*AP_max;
            cout << "steps: " << res.steps << endl;
            
            res.d_nu = 1/mes.T;
            
            cout << "d_nu: " << res.d_nu << endl;
            
            res.p_range.resize(res.steps);
            for (int i=0; i<res.steps; ++i)
                res.p_range[i] = i*res.d_nu;

            cout << "do bayes" << endl;
            bayesian_estimate_measures(mesP, comP, resP);
            
            
            write_measures(fileOut, comP, mesP, resP);
        }
        
       
        
	return 0;
}