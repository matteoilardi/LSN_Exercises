/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}
// Default constructor, does not perform any action

Random :: ~Random(){}
// Default destructor, does not perform any action

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   // This function generates a random number from a Gaussian distribution with given mean and sigma
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   // This function generates a random number in the range [min, max]
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  // This function generates a random number in the range [0,1)
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

double Random :: Exp(double lambda){
   // This function generates a random number from an exponential distribution with given mean (lambda)

   double x = Rannyu();
   return -1./lambda*log(1-x);
}

double Random :: Lorentz(double mean, double gamma){
   // This function generates a random number from an exponential distribution with given mean and width (gamma)

   double x = Rannyu();
   return gamma*tan(M_PI*(x-0.5)) + mean;
}

double Random :: circle_quarter(double min, double max){
   // This function generates a random number from a distribution with the shape of a quarter of a circle centered in the min, with radius max-min

   double r = max-min;
   double k = 4.*M_PI/(r*r);
   double x, y;
   do{
      x = Rannyu()*r;
      y = Rannyu()*r*k;
   }while(y > k*sqrt(r*r - x*x));

   return x + min;
}

int Random :: randInt(int min, int max){
   // This function generates a random integer in range (min, max)
   int a = int(Rannyu((double)min, (double)max+1));
   return a;
}

bool Random :: randBool(){
   // This function generates a random bool
   double a = Rannyu(0., 2.);
   return (a<1)? true : false;
}


double Random :: rand_theta(){
   // This function samples theta of a uniform distribution on the sphere
   double y = Rannyu(0., 1.);
   double theta = acos(1 - 2*y);
   return theta;
}

double Random :: rand_phi(){
   // This function samples phi of a uniform distribution on the sphere
   double phi = Rannyu(0., 2*M_PI);
   return phi;
}


void Random :: SetRandom(int * s, int p1, int p2){
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
