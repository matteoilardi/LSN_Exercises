#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "random.h"
#include "SA_Hamiltonian.h"

using namespace std;

const int nsteps_I = 1000;
const int nsteps_B = 10000;
double delta_I = 1.3;
double delta_B = 1.;
const int nblks = 100;

int main(){

    SA_Hamiltonian SA(nsteps_I, nsteps_B, delta_I, delta_B, nblks);

    // set initial parameters
    double mu0 = 0.9, sigma0 = 0.3; 
    SA.SetParams(mu0, sigma0);

    // annealing schedule
    // SA.ReadAnnealingSchedule("schedule.in");

    double gamma = 0.95;
    double T_0 = 10.;
    double T_min = 0.0001;
    SA.UniformCooling(gamma, T_0, T_min, 100);

    // simulation
    while(SA.GeticoolingStep() <= SA.GetncoolingSteps()){
        SA.CoolingStep();
        SA.Measure();
    }

    SA.FinalMeasurement();

    return 0;
}
