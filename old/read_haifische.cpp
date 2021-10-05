#include <iostream>
#include <netcdfcpp.h>

#include "read_haifische.h"

void read_simulation(string fileNet, network *netP, simulation *simP)
{
	cout << "reading " << fileNet << "..." << endl;
        
	NcFile ParaInNetwork(fileNet.c_str(), NcFile::ReadOnly);
	
// get the membrane constants
	NcVar* tau_AP = ParaInNetwork.get_var("tau_A");
	NcVar* tau_NP = ParaInNetwork.get_var("tau_N");
// 	NcVar* tau_GP = ParaInNetwork.get_var("tau_G");
	NcVar* tau_MP = ParaInNetwork.get_var("tau_M");
	
	tau_AP -> get(&netP->paras.tau_A,1);
	tau_NP -> get(&netP->paras.tau_N,1);
// 	tau_GP -> get(&netP->paras.tau_G,1);
	tau_MP -> get(&netP->paras.tau_M,1);
	
// get the interaction parameters of populations
	NcVar* kappaP = ParaInNetwork.get_var("kappa");
	NcVar* etaP = ParaInNetwork.get_var("eta");
	NcVar* epsP = ParaInNetwork.get_var("eps");
// 	NcVar* nP = ParaInNetwork.get_var("n");
	
	kappaP -> get(&netP->paras.kappa,1);
	etaP -> get(&netP->paras.eta,1);
	epsP -> get(&netP->paras.eps,1);
// 	nP -> get(&netP->paras.n,1);
// get simulation arrays
	NcVar* rateWntP = ParaInNetwork.get_var("rateWnt");
	NcVar* alpha_0P = ParaInNetwork.get_var("alpha_0");
	NcVar* tau_GP = ParaInNetwork.get_var("tau_G");
	NcVar* nP = ParaInNetwork.get_var("n");
	
	NcDim* rateWntSzP = ParaInNetwork.get_dim("rateWnt_dim");
	NcDim* alpha_0SzP = ParaInNetwork.get_dim("alpha_0_dim");
	NcDim* tau_GSzP = ParaInNetwork.get_dim("tau_G_dim");
	NcDim* nSzP = ParaInNetwork.get_dim("n_dim");
	
	int rateWntSz = rateWntSzP->size();
	int alpha_0Sz = alpha_0SzP->size();
	int tau_GSz = tau_GSzP->size();
	int nSz = nSzP->size();
	
	simP -> rateWnt.resize(rateWntSz);
	simP -> alpha_0.resize(alpha_0Sz);
	simP -> tau_G.resize(tau_GSz);
	simP -> n.resize(nSz);
	
	rateWntP -> get(&simP->rateWnt.front(),rateWntSz);
	alpha_0P -> get(&simP->alpha_0.front(),alpha_0Sz);
	tau_GP -> get(&simP->tau_G.front(),tau_GSz);
	nP -> get(&simP->n.front(),nSz);
}

void write_results(string fileOut, simulation sim, double *gammaP, double *chiP)
{
	cout << "writing " << fileOut << "..." << endl;
	
	NcFile writeResults(fileOut.c_str(), NcFile::Replace);
	
// 	cout << "size : " << gammaP.size() << endl;
// 	cout << "size: " << gammaP[0].size() << endl;
	
	NcDim* rateWnt_dim = writeResults.add_dim("rateWntSz", sim.rateWnt.size());
	NcDim* alpha_0_dim = writeResults.add_dim("alpha_0Sz", sim.alpha_0.size());
	NcDim* tau_G_dim = writeResults.add_dim("tau_GSz", sim.tau_G.size());
	NcDim* n_dim = writeResults.add_dim("nSz", sim.n.size());
	NcDim* pop_dim = writeResults.add_dim("populations", 2);
	
	NcVar *gamma = writeResults.add_var("gamma", ncDouble, n_dim, alpha_0_dim, tau_G_dim, rateWnt_dim, pop_dim);
	NcVar *chi = writeResults.add_var("chi", ncDouble, n_dim, alpha_0_dim, tau_G_dim, rateWnt_dim, pop_dim);

	gamma->put(gammaP, sim.n.size(), sim.alpha_0.size(), sim.tau_G.size(), sim.rateWnt.size(), 2);
	chi->put(chiP, sim.n.size(), sim.alpha_0.size(), sim.tau_G.size(), sim.rateWnt.size(), 2);
// 	
// 	chi->put(resP->chi, resP->gamma.size(), resP->gamma[0].size(), 2);
}