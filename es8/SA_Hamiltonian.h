#ifndef __SA_Hamiltonian__
#define __SA_Hamiltonian__

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "random.h"

using namespace std;

class SA_Hamiltonian{

    public:
    SA_Hamiltonian(int nsteps_I_in, int nsteps_B_in, double delta_I_in, double delta_B_in, int nblks_in){   
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
            while ( !input.eof() ){
                input >> property;
                if(property == "RANDOMSEED"){
                    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                    rnd.SetRandom(seed,p1,p2);
                }
            }
            input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

        icoolingStep = 1;
        
        nsteps_I = nsteps_I_in;
        nsteps_B = nsteps_B_in;
        delta_I = delta_I_in;
        delta_B = delta_B_in;
        nblks = nblks_in;

        nbins = 300;
        waveFunction = new double[nbins];
        for(int ibin = 0; ibin < nbins; ibin++){
            waveFunction[ibin] = 0;
        }
        start = -3.;
        end = 3.;
        npoints = 0;
    };
    ~SA_Hamiltonian(){delete[] waveFunction;};

    void SetAnnealingSchedule(vector<double> temps, vector<int> coolingSteps){
        annealingSchedule_steps = coolingSteps;
        annealingSchedule_temps = temps;

        ncoolingSteps = temps.size();
    };

    void ReadAnnealingSchedule(string fileName){
        double T;
        int steps;
        ifstream In;
        In.open(fileName);
        
        In >> T >> steps;
        while(!In.eof()){
            annealingSchedule_temps.push_back(T);
            annealingSchedule_steps.push_back(steps);
            In >> T >> steps;
        }

        ncoolingSteps = annealingSchedule_temps.size();
        In.close();
    }

    void UniformCooling(double coolingRate, double T_0, double T_min, double steps){
        double T = T_0;
        while(T > T_min){
            annealingSchedule_temps.push_back(T);
            annealingSchedule_steps.push_back(steps);
            T = T*coolingRate;
        }

        ncoolingSteps = annealingSchedule_temps.size();
    }

    void SetParams(double mu0, double sigma0){
        mu = mu0;
        sigma = sigma0;
    };

    int GeticoolingStep() const{
        return icoolingStep;
    };

    int GetncoolingSteps() const{
        return ncoolingSteps;
    };


    void CoolingStep(){
        if(icoolingStep > ncoolingSteps){
            cerr << "Annealing schedule completed!" << endl;
            return;
        } 
    
        double mu_proposed, sigma_proposed;
        double nsteps = annealingSchedule_steps[icoolingStep-1];
        double T = annealingSchedule_temps[icoolingStep-1];
        double attempted_mu = 0., accepted_mu = 0., attempted_sigma = 0., accepted_sigma = 0.;

    
        for(int istep = 0; istep < nsteps; istep++){
            mu_proposed = mu + rnd.Rannyu(-delta_B, delta_B);
            attempted_mu += 1.;
            if(rnd.Rannyu() < Boltzmann(T, mu_proposed, sigma)/Boltzmann(T, mu, sigma)){ //conditional probability for sigma
                mu = mu_proposed;
                accepted_mu += 1.;
            }

            sigma_proposed = sigma + rnd.Rannyu(-delta_B, delta_B);
            attempted_sigma += 1.;
            if(rnd.Rannyu() < Boltzmann(T, mu, sigma_proposed)/Boltzmann(T, mu, sigma)){ //conditional probability for sigma
                sigma = sigma_proposed;
                accepted_sigma += 1.;
            }
        }

        cout << endl << "---------------------------------------------------------" << endl << endl;
        cout << "Cooling step #" << icoolingStep << endl;
        cout << "Acceptance rate for mu and sigma: " << accepted_mu/attempted_mu << ", " << accepted_sigma/attempted_sigma << endl;

        if(accepted_mu/attempted_mu < 0.35 or accepted_sigma/attempted_sigma < 0.35){
            delta_B = delta_B/5.; //keep the acceptance rate high enough for the Metropolis algorithm to work

            cout << "Upgraded Metropolis step = " << delta_B << endl;
        }
        icoolingStep++;
    };

    void Measure(){
        double sum = 0.;
        double sum2 = 0;
        double energy;

        for(int iblk = 1; iblk <= nblks; iblk++){
            energy = trial_energy(mu, sigma, false);

            sum += energy;
            sum2 += energy*energy;
        }
        err_H = sum2/nblks - pow(sum/nblks, 2);
        err_H = sqrt(err_H);
        H = sum/nblks;

        ofstream Energy;
        Energy.open("energy.out",ios::app);
        Energy << icoolingStep-1 << '\t' << H << '\t' << err_H << endl;
        Energy.close();

        ofstream Params;
        Params.open("params.out", ios::app);
        Params << icoolingStep-1 << '\t' << mu << '\t' << sigma << endl;
        Params.close();
    };

    void FinalMeasurement(){
    
        double T = annealingSchedule_temps[icoolingStep-2];
        cout << endl;
        cout << "Final temperature = " << T << endl;

        // final values of the parameters mu and sigma and uncertainties
        cout << setprecision(3) << "Parameters: mu = " << mu << ", sigma = " << sigma << endl;


        // block average of the trial state energy
        cout << endl;
        cout << "Meausring the energy of the optimal ground state..." << endl;

        ofstream EnergyFinal;
        EnergyFinal.open("energy_fin.out",ios::app);
        double sum_H = 0., sum2_H = 0.;
        for(int iblk = 1; iblk <= nblks; iblk++){
            H = trial_energy(mu, sigma, true);
            sum_H += H;
            sum2_H += H*H;

            err_H = sum2_H/iblk - pow(sum_H/iblk, 2);
            err_H = sqrt(err_H);
            H = sum_H/iblk;
            EnergyFinal << iblk << '\t' << H << '\t' << err_H/sqrt(iblk-1) << endl;
        }

        EnergyFinal.close();

        // psi2
        cout << "Filling the histogram for the square modulus of the wave function..." << endl;

        ofstream Psi2;
        double binSize, xbin, ybin;
        binSize = (end-start)/nbins;

        Psi2.open("psi2.out", ios::app);
        for(int ibin = 0; ibin < nbins; ibin++){
            xbin = -3. + binSize/2. + binSize*ibin;
            ybin = waveFunction[ibin]/npoints/binSize;
            Psi2 << xbin << '\t' << ybin << endl;
        }

        Psi2.close();

        cout << "Done. Data in files energy.out, energy_fin.out, params.out, psi2.out." << endl;
        cout << endl;
    }


    double Boltzmann(double T, double mu, double sigma){
        double energy = trial_energy(mu, sigma, false);
        double beta = 1./T;

        return exp(-1.*beta*energy);
    };

    double trial_energy(double mu, double sigma, bool wf){
        double I = 0;
        double x = rnd.Rannyu(-delta_I, delta_I);
        double y;

        for(int istep = 0; istep < nsteps_I; istep++){
            y = x + rnd.Rannyu(-delta_I, delta_I);
        
            if(rnd.Rannyu() < psi2(y, mu, sigma)/psi2(x, mu, sigma)){
                x = y;
            }


            I += E_loc(x, mu, sigma)/(double)nsteps_I;

            if(wf){
                int bin = nbins * (x+3.)/(end-start);
                if(0 <= bin and bin <= nbins){
                    waveFunction[bin]++;
                    npoints++;
                }
            }
        }

        return I;
    }

    double E_loc(double x, double mu, double sigma){
        // hbar = 1, m  = 1
        double delta1 = pow(x - mu, 2);
        double delta2 = pow(x + mu, 2);
        double Tpsi = -1./(2*pow(sigma, 4)) * ( exp(-delta1/(2*sigma*sigma))*(delta1 - sigma*sigma) + exp(-delta2/(2*sigma*sigma))*(delta2 - sigma*sigma) );
        double psi = exp(-delta1/(2*sigma*sigma)) + exp(-delta2/(2*sigma*sigma));
        double V = pow(x, 4) - 5./2. * x*x;

        return Tpsi/psi + V;
    }

    double psi2(double x, double mu, double sigma){
        double delta1 = pow(x - mu, 2);
        double delta2 = pow(x + mu, 2);
        double psi = exp(-delta1/(2*sigma*sigma)) + exp(-delta2/(2*sigma*sigma));

        return psi*psi;
    }

    double V(double x){
        return pow(x, 4) - 5./2.*pow(x, 2);
    };

    
    private:
    Random rnd;

    vector<int> annealingSchedule_steps;
    vector<double> annealingSchedule_temps;
    int icoolingStep, ncoolingSteps;

    
    int nblks; // data blocking for the energy of the trial state at each cooling step, for mu and sigma at the end

    double delta_I; // Metropolis step length for integral evaluation
    int nsteps_I; // #steps for integral evaluation

    double delta_B; // Metropolis step length for Boltzmann distribution (both for mu and sigma: this is reasonable if we assume they are the same order of magnitude)
    int nsteps_B; // block length in the final calculation of mu and sigma

    double H, err_H;  // energy
    double mu, sigma; // optimization parameters

    double* waveFunction; // wave function histogram
    int nbins, npoints;
    double start, end; 
};

#endif // __SA_Hamiltonian__
