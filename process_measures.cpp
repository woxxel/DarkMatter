#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdio.h>

#include "bayes_estimate.cpp"

using namespace std;

int main()
{
        if ((argc < 4) or (argc > 4))
        {
                cout << "Please specify (only) the input and output file name!" << endl;
                return 1;
        }
        
        string fileNet = argv[1];
        string fileSim = argv[2];
        string fileOut = argv[3];

        network net, *netP = &net;
        simulation sim, *simP = &sim;

        results res, *resP = &res;
        
        read_measures(fileNet,netP,fileSim,simP);
        
        
}