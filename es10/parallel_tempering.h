#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "random.h"
#include "mpi.h"

using namespace std;

// loaded dice parameter
#define P1 5.

class Map{
    public:
    Map(string fileName);
    ~Map(){};

    int get_ncities(){return m_ncities;};
    double get_x(int icity){return m_x[icity];};
    double get_y(int icity){return m_y[icity];};
    double distance1(int i, int j){
        return dist_matrix[i][j];
    };
    double distance2(int i, int j){
        return pow(dist_matrix[i][j], 2);
    };

    private:
    int m_ncities;
    vector<double> m_x;
    vector<double> m_y;
    vector<vector<double>> dist_matrix;
};



class Path{
    public:
    Path(Map* pmap, double T, int rank, int size);
    ~Path(){};

    double get_length(int n) const{
        if(n == 1) return m_length1;
        if(n == 2) return m_length2;

        return 0;
    };

    int get_city(int icity){return m_path[icity];};
    double get_temperature(){return m_T;};
    vector<int> get_path(){return m_path;};

    void set_path(int* new_path){
        for(int i = 0; i < m_ncities; i++) m_path[i] = new_path[i];
    }

    int* to_array(){
        int* array = new int[m_ncities];

        for(int i = 0; i < m_ncities; i++) array[i] = get_city(i);
        return array;
    };

    
    // checks if the path is legitimate, i. e. every city appears once and the city #0 is the first one 
    bool check();
    // stores the path's length in the appropriate data member
    void measure();
    // prints to a file the coordinates of the cities in the path
    void print(string fileName);

    // Metropolis moves: genetic mutation operators
    void permutation();
    void shift();
    void block_swap();
    void inversion();

    // Metropolis move involving the swap with the first path at an higher temperature
    void temperature_swap_intra(Path& path_up);
    // Metropolis move involving path swaps between different processes
    void temperature_swap_inter(bool start);
    // Monte Carlo step
    void mutate();

    // computes the probability of accepting the Monte Carlo move
    double p_Boltzmann(Path trial_path);
    // computes the probability of swapping two paths at different temperatures
    double p_swap_T_Boltzmann(double trial_length, double trial_T);
    
    

    // periodic boundary conditions
    int pbc(int icity);
    // tells whether city #number appears in the path
    bool is_in(int number);
    // draws a random int - small numbers are more likely
    int loaded_dice(double p = P1);

    bool operator<(const Path& other) const{
        return get_length(1) < other.get_length(1);
    }
    
    private:
    Random m_rnd;
    Map* m_pmap;
    
    vector<int> m_path;
    int m_ncities;
    double m_length1, m_length2;

    double m_T;
    int m_rank, m_size;
};

// Metropolis moves involving temperature swaps within the same process
void temperature_swaps_intra(vector<Path> paths);

// Reads from a file the temperatures used for parallel tempering
vector<double> readTemperatures(string fileName);
