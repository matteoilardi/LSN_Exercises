#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

#define M 1000000 
#define N 100 // # blocks
#define n_intervals 100

using namespace std;
// generates one random value of the call option at t = 0 directly sampling the price
double call_option_price_c(Random* prnd, double S_0, double T, double K, double r, double sigma);
// generates one random value of the call option at t = 0 dividing time evolution in n_intervals steps
double call_option_price_d(Random* prnd, double S_0, double T, double K, double r, double sigma);
// generates one random value of the put option at t = 0 directly sampling the price
double put_option_price_c(Random* prnd, double S_0, double T, double K, double r, double sigma);
// generates one random value of the put option at t = 0 dividing time evolution in n_intervals steps
double put_option_price_d(Random* prnd, double S_0, double T, double K, double r, double sigma);
double max(double x, double y);

double mean(vector<double>& v, unsigned int j);
double stddev(vector<double>& v, unsigned int j);
void print_progressive_mean_stdev(ofstream* pOut, vector<double> v);


int main(int argc, char *argv[]){

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


    const double S_0 = 100.; // initial asset price
    const double T = 1.; // delivery time
    const double K = 100.; // strike price
    const double r = 0.1; // risk-free interest rate
    const double sigma = 0.25; // volatility

    const int L = M/N; // # points per block
    double block_mean_cc;
    double block_mean_cd;
    double block_mean_pc;
    double block_mean_pd;

    vector<double> call_option_c;
    vector<double> call_option_d;
    vector<double> put_option_c;
    vector<double> put_option_d;


    for(int i = 0; i < N; i++){
        block_mean_cc = 0;
        block_mean_cd = 0;
        block_mean_pc = 0;
        block_mean_pd = 0;
        
        for(int k = 0; k < L; k++){
            block_mean_cc += call_option_price_c(&rnd, S_0, T, K, r, sigma)*1./L;
            block_mean_cd += call_option_price_d(&rnd, S_0, T, K, r, sigma)*1./L;
            block_mean_pc += put_option_price_c(&rnd, S_0, T, K, r, sigma)*1./L;
            block_mean_pd += put_option_price_d(&rnd, S_0, T, K, r, sigma)*1./L;
        }

        call_option_c.push_back(block_mean_cc);
        call_option_d.push_back(block_mean_cd);
        put_option_c.push_back(block_mean_pc);
        put_option_d.push_back(block_mean_pd);
    }


    ofstream Out1("call_option_c.out");
    ofstream Out2("call_option_d.out");
    ofstream Out3("put_option_c.out");
    ofstream Out4("put_option_d.out");

    print_progressive_mean_stdev(&Out1, call_option_c);
    print_progressive_mean_stdev(&Out2, call_option_d);
    print_progressive_mean_stdev(&Out3, put_option_c);
    print_progressive_mean_stdev(&Out4, put_option_d);


    cout << endl;
    cout << "Data (# blocks, estimation, error) in files:" << endl;
    cout << "call_option_c.out" << endl;
    cout << "call_option_d.out" << endl;
    cout << "put_option_c.out" << endl;
    cout << "put_option_d.out" << endl;
    cout << endl;
    

    Out1.close();
    Out2.close();
    Out3.close();
    Out4.close();
    rnd.SaveSeed();
    return 0;
}


double call_option_price_c(Random* prnd, double S_0, double T, double K, double r, double sigma){
    double W = prnd->Gauss(0., T);
    double S_T = S_0 * exp((r - 0.5*sigma*sigma)*T + sigma*W);

    return max(0., S_T-K)*exp(-1.*r*T);
}

double call_option_price_d(Random* prnd, double S_0, double T, double K, double r, double sigma){
    double Z;
    double growth = 1.;
    double deltat = T/n_intervals;
    double S_T;

    for(int i = 0; i < n_intervals; i++){
        Z = prnd->Gauss(0., 1.);
        growth = growth * exp((r - 0.5*sigma*sigma)*deltat + sigma*Z*sqrt(deltat));
    }

    S_T = S_0 * growth;
    return max(0., S_T-K)*exp(-1.*r*T);
}

double put_option_price_c(Random* prnd, double S_0, double T, double K, double r, double sigma){
    double W = prnd->Gauss(0., T);
    double S_T = S_0 * exp((r - 0.5*sigma*sigma)*T + sigma*W);

    return max(0., K-S_T)*exp(-1.*r*T);
}

double put_option_price_d(Random* prnd, double S_0, double T, double K, double r, double sigma){
    double Z;
    double growth = 1.;
    double deltat = T/n_intervals;
    double S_T;

    for(int i = 0; i < n_intervals; i++){
        Z = prnd->Gauss(0., 1.);
        growth = growth * exp((r - 0.5*sigma*sigma)*deltat + sigma*Z*sqrt(deltat));
    }

    S_T = S_0 * growth;
    return max(0., K-S_T)*exp(-1.*r*T);
}

double max(double x, double y){
    return (x>y)? x:y;
}


void print_progressive_mean_stdev(ofstream* pOut, vector<double> v){
   for(int j = 1; j < v.size(); j++){
      *pOut << j+1 << '\t' << left << setw(12) << mean(v, j) << setw(12) << stddev(v, j) << endl;
   }
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
