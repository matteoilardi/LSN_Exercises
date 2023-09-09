#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "MonteCarlo.h"

using namespace std;

distribution::distribution(double min, double max){
    m_min = min;
    m_max = max;

    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    }else cerr << "PROBLEM: Unable to open Primes" << endl;
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

    rnd.SaveSeed(); //non ha senso!
}

double MC(distribution* p, double (*f)(double), int n_random_block){
    double x;
    double sum = 0;

    for(int i = 0; i < n_random_block; i++){
        x = p->Rand();
        sum += f(x)/(p->Eval(x)*n_random_block);
    }

    return sum;
}
