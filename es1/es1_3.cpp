#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"

#define M 10000 // throws
#define N 100 // blocks

using namespace std;

// generates a random number from a cosine distribution in (0, 1) without using the exact value of pi
double rand_cos(Random* prnd);

void print_progressive_mean_stdev(ofstream* pOut, vector<double> v);
double mean(vector<double>& v, unsigned int j);
double stddev(vector<double>& v, unsigned int j);


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


   const double l = 1.; //stick length
   const double d = 3.; //line spacing
   vector<double> pi;
   const int L = M/N; // throws per block

   int count; // ++ if stick intersects the line
   double cos_theta, h;
   double p;

   for(int i = 0; i < N; i++){
      count = 0;

      for(int j = 0; j < L; j++){
         h = rnd.Rannyu(0., d);
         cos_theta = rand_cos(&rnd);
         count += (h + l*cos_theta) > d;
      }

      p = static_cast<double>(count)/L;
      pi.push_back(2*l/(p*d));
   }

   // Progressive mean
   ofstream Out("Buffon.out");
   print_progressive_mean_stdev(&Out, pi);
   Out.close();

   cout << endl;
   cout << "Results in Buffon.out. # blocks, mean, standard deviation";
   cout << endl;

   rnd.SaveSeed();
   return 0;
}


double rand_cos(Random* prnd){
    double x, y;
    do{
        x = prnd->Rannyu();
        y = prnd->Rannyu();
    }while(sqrt(x*x + y*y) > 1);

    return x/sqrt(x*x + y*y);
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
