using namespace std;

void read_model(string fileModel, model *modP);
void read_computation(string fileComp, computation *comP);
void read_measures(string fileMeasures, measures *mesP);

// void write_results(string fileOut, computation sim, double *gammaP, double *chiP);
void write_theory(string fileOut, computation *comP, results *resP);
void write_theory(string fileOut, computation *comP, measures *mesP, results *resP);