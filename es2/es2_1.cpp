#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "MonteCarlo.h"

#define M 10000
#define N 100  // maximum number of blocks

using namespace std;

double f(double x);
double mean(vector<double>& v, unsigned int j);
double stddev(vector<double>& v, unsigned int j);
void print_progressive_mean_stdev(ofstream* pOut, vector<double> v);


int main (int argc, char *argv[]){
   const int L = M/N; // random numbers per block

   // Uniform distribution
   uniform unf(0., 1.);
   vector<double> I_unif;
   double I;

   for(int i = 0; i < N; i++){
      I = MC(&unf, f, L);
      I_unif.push_back(I);
   }

   // Circle-quarter-like distribution
   circle_quarter cq(0., 1.);
   vector<double> I_circle;

   for(int i = 0; i < N; i++){
      I = MC(&cq, f, L);
      I_circle.push_back(I);
   }

   // Linear distribution
   linear lin(0., 1., 0.1);
   vector<double> I_linear;

   for(int i = 0; i < N; i++){
      I = MC(&lin, f, L);
      I_linear.push_back(I);
   }


   // Progressive mean
   ofstream Out1("int_uniform.out");
   ofstream Out2("int_circle.out");
   ofstream Out3("int_linear.out");


   print_progressive_mean_stdev(&Out1, I_unif);
   print_progressive_mean_stdev(&Out2, I_circle);
   print_progressive_mean_stdev(&Out3, I_linear);

   Out1.close();
   Out2.close();
   Out3.close();

   cout << endl;
   cout << "Data in int_uniform.out, int_circle.out, int_linear.out: # blocks, mean, standard deviation.";
   cout << endl;

   return 0;
}




double f(double x){
   double y = 0.5*M_PI*cos(0.5*M_PI*x);
   return y;
};

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
