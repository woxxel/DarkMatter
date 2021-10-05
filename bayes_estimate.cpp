#include <vector>
#include <math.h>

void bayesian_estimate_theory(model *modP, computation *comP, results *resP)
{
        // wanna calculate: P(lambda|k) = P(k|lambda)*p(lambda)/p(k)
        // if no prior beliefs exist, p(lambda) = p(k) -> P(lambda|k) = P(k|lambda)
        // however, if we take p(lambda) to be firing rate distribution from theory, we obtain the following:
        
//     P_k_lambda  = poisson_distr(k,nu*T);
//     p_lambda    = modP->distribution_exact(nu,0);
//     p_k = int(P_k_lambda*p_lambda) d lambda
        
        resP->p_bayes_est[comP->k][comP->j].resize(resP->steps,0);
        for (int n=0; n<comP->N; ++n)
        {
//                 cout << "APs for neuron " << n << ": " << k << endl;
                double p_k = 0;
                
                for (int i=1; i < resP->steps; ++i)    // "integration" over lambda    //should also work with integration (faster?)
                        p_k += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior);
                
                resP->p_bayes_est[comP->k][comP->j][0] = 0;
                for (int i=1; i < resP->steps; ++i)    // construction of posterior
                        resP->p_bayes_est[comP->k][comP->j][i] += poisson_distr(resP->N_AP[comP->k][comP->j][n],resP->factor*i)*bayes_est_prior(i*resP->d_nu,modP,comP->prior)/p_k / (resP->d_nu*comP->N);
        }
}


void bayesian_estimate_measures(measures *mesP, computation *comP, results *resP)
{
        // wanna calculate: P(lambda|k) = P(k|lambda)*p(lambda)/p(k)
        // if no prior beliefs exist, p(lambda) = p(k) -> P(lambda|k) = P(k|lambda)
        // however, if we take p(lambda) to be firing rate distribution from theory, we obtain the following:
        
//     P_k_lambda  = poisson_distr(k,nu*T);
//     p_lambda    = modP->distribution_exact(nu,0);
//     p_k = int(P_k_lambda*p_lambda) d lambda
        
        resP->p_bayes_est_measures.resize(resP->steps,0);
        for (int n=0; n<mesP->N; ++n)
        {
//                 cout << "APs for neuron " << n << ": " << k << endl;
                double p_k = 0;
                
                for (int i=1; i < resP->steps; ++i)    // "integration" over lambda    //should also work with integration (faster?)
                        p_k += poisson_distr(mesP->N_AP[n],i)*bayes_est_prior_measures(i*resP->d_nu,comP->prior);
                
                resP->p_bayes_est_measures[0] = 0;
                for (int i=1; i < resP->steps; ++i)    // construction of posterior
                        resP->p_bayes_est_measures[i] += poisson_distr(mesP->N_AP[n],i)*bayes_est_prior_measures(i*resP->d_nu,comP->prior)/p_k / (resP->d_nu*mesP->N);
        }
}



double bayes_est_prior(double nu, model *modP, string prior)
{
        if (not prior.compare("mean_field"))
            return modP->distribution_exact(nu,0);
        else if (not prior.compare("plain"))
            return 1;
        else
        {
            cout << "Specify a prior!" << endl;
            exit(0);
        }
}


double bayes_est_prior_measures(double nu, string prior)
{
//         if (not prior.compare("mean_field"))
//             return modP->distribution_exact(nu,0);
        if (not prior.compare("plain"))
            return 1;
        else
        {
            cout << "Specify a prior!" << endl;
            exit(0);
        }
}