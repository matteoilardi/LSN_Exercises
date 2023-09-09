#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

#define N 10000

using namespace std;

int main(int argc, char *argv[]){

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;


    // four vectors for each distribution, storing the N = 10000 realizations of S_n, n = 1, 2, 10, 100
    const vector<int> n_block = {1, 2, 10, 100};
    double x;
    double sum;
    vector<double> Uniform[n_block.size()];
    vector<double> Exponential[n_block.size()];
    vector<double> Cauchy_Lorentz[n_block.size()];


    // uniform distribution
    for(int i = 0; i < n_block.size(); i++){
        for(int j = 0; j < N; j++){
            sum = 0;

            for(int k = 0; k < n_block[i]; k++){
                x = rnd.Rannyu();
                sum += x/n_block[i];
            }

            Uniform[i].push_back(sum);
        }
    }

    // exponential distibution
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < N; j++){
            sum = 0;

            for(int k = 0; k < n_block[i]; k++){
                x = rnd.Exp(1);
                sum += x/n_block[i];
            }

            Exponential[i].push_back(sum);
        }
    }

    // lorentzian distribution
   for(int i = 0; i < 4; i++){
        for(int j = 0; j < N; j++){
            sum = 0;

            for(int k = 0; k < n_block[i]; k++){
                x = rnd.Lorentz(0, 1);
                sum += x/n_block[i];
            }

            Cauchy_Lorentz[i].push_back(sum);
        }
    }

    cout << endl;
    cout << "Mean values of n = 1, 2, 10, 100 random numbers in uniform.out, exponential.out, lorentzian.out." << endl;
    cout << endl;

    ofstream Out1("uniform.out");
    ofstream Out2("exponential.out");
    ofstream Out3("lorentzian.out");

    for(int i = 0; i < N; i++){
        Out1 << left << setw(15) << Uniform[0][i] << setw(15) << Uniform[1][i] << setw(15) << Uniform[2][i] << setw(15) << Uniform[3][i] << endl;
    }

    for(int i = 0; i < N; i++){
        Out2 << left << setw(15) << Exponential[0][i] << setw(15) << Exponential[1][i] << setw(15) << Exponential[2][i] << setw(15) << Exponential[3][i] << endl;
    }

    for(int i = 0; i < N; i++){
        Out3 << left << setw(15) << Cauchy_Lorentz[0][i] << setw(15) << Cauchy_Lorentz[1][i] << setw(15) << Cauchy_Lorentz[2][i] << setw(15) << Cauchy_Lorentz[3][i] << endl;
    }


    Out1.close();
    Out2.close();
    Out3.close();
    rnd.SaveSeed();
    return 0;
}
