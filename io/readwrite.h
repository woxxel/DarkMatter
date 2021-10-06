using namespace std;

void get_from_ncid(int ncid, const char *varName, void *varP);
void get_from_ncid(int ncid, const char *varName, size_t *varP);

void write_prep_paras(int ncid, int dimSz, int *dimID, simulation *simP);
void write_paras(int ncid, int *dimID, simulation *simP);

// void write_to_ncid(int ncid, const char *varName, nc_type varType, int dimSz, const int *dimids, double *varData);
// void write_to_ncid(int ncid, const char *varName, nc_type varType, int dimSz, const int *dimids, double *varData, hyperslab hslab);

void read_model(string fileModel, model *modP);
void read_simulation(string fileSim, simulation *simP);
void read_computation(string fileComp, computation *comP);
void read_measures(string fileMeasures, measures *mesP);

// void write_results(string fileOut, computation sim, double *gammaP, double *chiP);
void write_theory(string fileOut, results *resP);
// void write_theory(string fileOut, computation *comP, measures *mesP, results *resP);
void write_sharks(string fileOut, simulation *simP, model *modP, results *resP, double *gammaP, double *chiP, double *regionsP);

void write_stats(string fileOut, model *modP, simulation *simP, results *resP);

// void write_stats_mode_3(NcFile writeResults, results *resP, int p, string addStr, int dim1, int dim2);

string strApp(string str1, string str2);
