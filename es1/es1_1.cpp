#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

#define M 10000
#define N 100  // total number of blocks
#define chi2_intervals 100

using namespace std;

double mean(vector<double>& v, unsigned int j);
double stddev(vector<double>& v, unsigned int j);


int main (int argc, char *argv[]){

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



   const int L = M/N; // random numbers per block
   double x = 0;
   int y;
   double sum, sum2;
   vector<double> A; // will store the difference between the average of x and the expected value 0.5 for each block
   vector<double> averages;
   vector<double> errors;

   vector<double> A2; // will store the difference between the average of sigma^2 and the expected value 1/12 squared for each block
   vector<double> averages2;
   vector<double> errors2;

   vector<int> chi_obs(chi2_intervals, 0);
   vector<double> chi2(N, 0.);
   
   ofstream Out1("average.out");
   ofstream Out2("variance.out");
   ofstream Out3("chi2.out");


   for(int i = 0; i < N; i++){
      sum = 0;
      sum2 = 0;

      for(int k = 0; k<L; k++){
         x = rnd.Rannyu();
         sum += x/L; 
         sum2 += pow(x-0.5, 2)/L;

         y = floor((x/1.)*chi2_intervals);
         chi_obs[y%chi2_intervals]++;
      }

      A.push_back(sum);
      A2.push_back(sum2);

      double chi_exp = (double)L*(i+1)/chi2_intervals;
      for(int j = 0; j < chi2_intervals; j++){
         chi2[i] += pow( static_cast<double>(chi_obs[j]) - chi_exp, 2)/chi_exp;
      }
   }

   for(int j = 1; j < N; j++){
      averages.push_back(mean(A, j));
      errors.push_back(stddev(A, j));
      Out1 << (j+1)*L << '\t' << averages[j-1] << '\t' << errors[j-1] << endl;

      averages2.push_back(mean(A2, j));
      errors2.push_back(stddev(A2, j));
      Out2 << (j+1)*L << '\t' << averages2[j-1] << '\t' << errors2[j-1] << endl;

      Out3 << (j+1)*L << '\t' << chi2[j] << endl;
   }

   cout << endl;
   cout << "Data in average.out, variance.out and chi2.out." << endl;
   cout << endl;

   Out1.close();
   Out2.close();
   Out3.close();
   rnd.SaveSeed();
   return 0;
}



// returns the average of the first (j+1) elements of the vector v
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

// returns the standard deviation of the first (j+1) elements of the vector v
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
