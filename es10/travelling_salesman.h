#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "random.h"

using namespace std;

// genetic algorithm parameters
#define P1 5.
#define P2 0.08
#define P3 0.1
#define P4 0.1
#define P5 0.1
#define P6 5.

class Map{
    public:
    Map(string fileName, Random* pr);
    Map(int ncities, int type, Random* pr);
    ~Map(){};

    int get_ncities(){return m_ncities;};
    double get_x(int icity){return m_x[icity];};
    double get_y(int icity){return m_y[icity];};
    string get_type(){return typeName;};
    double distance1(int i, int j){
        return dist_matrix[i][j];
    };
    double distance2(int i, int j){
        return pow(dist_matrix[i][j], 2);
    };

    private:
    Random* m_pr;

    string typeName;
    int m_ncities;
    vector<double> m_x;
    vector<double> m_y;
    vector<vector<double>> dist_matrix;
};



class Path{
    public:
    Path(Map* pmap, Random* pr);
    ~Path(){};

    double get_length(int n) const{
        if(n == 1) return m_length1;
        if(n == 2) return m_length2;

        return 0;
    };
    vector<int>& get_path(){return m_path;};
    int get_city(int icity){return m_path[icity];};
    
    // checks if the path is legitimate, i. e. every city appears once and the city #0 is the first one 
    bool check();
    // stores the path's length in the appropriate data member
    void measure();
    // prints to a file the coordinates of the cities in the path
    void print(string fileName);
    // returns the name of the map type
    string get_map_type();

    // genetic mutation operators
    void permutation();
    void shift();
    void block_swap();
    void inversion();

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
    Random* m_pr;
    Map* m_pmap;
    
    vector<int> m_path;
    int m_ncities;
    double m_length1, m_length2;
};



class Population{
    public:
    Population(Map& map, int npaths, Random* pr);
    ~Population(){};

    // introduces random mutations in the population
    void mutate();
    // generates npaths new individuals and preserves the npaths most fit, choosing both from parents and the offspring
    void newgen();
    // adds to the population two children generated from two given parents
    void crossover(int ipath, int jpath);

    // sorts paths by length (l1)
    void sort_paths();
    // returns the average length (l1) of the first half of the population
    double half_av_l1();
    // print the average length (l1) of the first half of the population
    void print_half_av_l1(int igen);
    // returns a pointer to the shortest (l1) path in the population
    Path get_best_path_l1();


    private:
    Random* m_pr;

    int m_npaths;
    int m_ncities;
    vector<Path> m_pop;
};