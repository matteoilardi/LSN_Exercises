#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "RandomWalk.h"

using namespace std;

#define M 10000 // number of random walks
#define N 100 // number of blocks
#define n_steps 100

double mean(vector<double>& v, unsigned int j);
double stddev(vector<double>& v, unsigned int j);
void print_progressive_mean_stdev(ofstream* pOut, vector<double> v);


int main(int argc, char *argv[]){

    // RNG setup
    Random rnd;
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

    ofstream Out1("RW_lattice.out");
    ofstream Out2("RW_continuum.out");


    int L = M/N; // random walks per block

    // random walk initialization
    vector<RW_lattice> RW_l;
    vector<RW_continuum> RW_c;


    RW_lattice rwl(1., &rnd);
    RW_continuum rwc(1., &rnd);
    for(int i = 0; i < M; i++){
        RW_l.push_back(rwl);
        RW_c.push_back(rwc);
    }
    

    for(int i = 0; i < n_steps; i++){
        vector<double> vl(N, 0.);
        vector<double> vc(N, 0.);

        // calculates the square root of <r^2> in each block 
        for(int j = 0; j < M; j++){
            RW_l[j].step();
            RW_c[j].step();

            vl[j/L] += RW_l[j].r2()/L;
            vc[j/L] += RW_c[j].r2()/L;

            if(! ((j+1)%L)){
                vl[j/L] = sqrt(vl[j/L]);
                vc[j/L] = sqrt(vc[j/L]);
            }
        }

        // average and standard deviation among blocks, at each step
        Out1 << i+1 << '\t' << left << setw(12) << mean(vl, N-1) << setw(12) << stddev(vl, N-1) << endl;
        Out2 << i+1 << '\t' << left << setw(12) << mean(vc, N-1) << setw(12) << stddev(vc, N-1) << endl;

    }

    cout << endl;
    cout << "Results in RW_lattice.out and RW_continuum.out. # steps, mean, standard deviation.";
    cout << endl;

    Out1.close();
    Out2.close();
    rnd.SaveSeed();
    return 0;
}


double mean(vector<double>& v, unsigned int j){
    unsigned int n = j+1;
	double sum = 0.;

    if(v.begin()+j > v.end()){
        cerr << "Vector length exceeded!" << endl;
        return 0.;
    }

    for(auto i = v.begin(); i <= v.begin()+j; i++){
        sum += (*i)/n;
    }

	return sum;
};

double stddev(vector<double>& v, unsigned int j){
    unsigned int n = j+1;
    double sum = 0.;
    double sum2 = 0.;

    if(v.begin()+j > v.end()){
        cerr << "Vector length exceeded!" << endl;
        return 0.;
    }

    for(auto i = v.begin(); i <= v.begin()+j; i++){
        sum += (*i)/n;
        sum2 += (*i)*(*i)/n;
    }

    return sqrt(sum2 - sum*sum)/sqrt(n-1);
};
