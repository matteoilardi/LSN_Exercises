#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "random.h"
#include "travelling_salesman.h"

using namespace std;



int main(){

    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if(Primes.is_open()){
        Primes >> p1 >> p2 ;
    }else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if(input.is_open()){
        while (!input.eof()){
            input >> property;
            if(property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }else cerr << "PROBLEM: Unable to open seed.in" << endl;



    const int ncities = 34;
    const int npaths = 1000;
    const int ngenerations = 1000;

    vector<double> average_l1;

    Map map(ncities, 0, &rnd);
    Population population(map, npaths, &rnd);

    for(int igen = 1; igen <= ngenerations; igen++){
        population.mutate();
        population.sort_paths();
        population.newgen();

        population.print_half_av_l1(igen);
    }

    Path best_path = population.get_best_path_l1();
    best_path.print("best_path");

    cout << endl;
    cout << "Length (l1) of the shortest path = " << best_path.get_length(1) << endl;
    cout << endl;

    return 0;
}
