#include <iostream>
#include <string>
#include <fstream>
#include "parallel_tempering.h"

using namespace std;

int main(int argc, char* argv[]){
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<double> temperatures = readTemperatures("temperatures.in");
    if(temperatures.size()%size != 0){
        if(rank == 0) cerr << endl << "The number of temperatures must be a multiple of the number of processes!" << endl;
        return 1;
    }
    int ntemps = temperatures.size()/size;

    Map American_capitals("American_capitals.in");
    vector<Path> paths;
    for(int itemp = 0; itemp < ntemps; itemp++){
        Path path(&American_capitals, temperatures[rank*ntemps + itemp], rank, size);
        paths.push_back(path);
    }

    const int nsteps = 10000; 

    for(int istep = 0; istep < nsteps; istep++){
        for(int itemp = 0; itemp < ntemps; itemp++){

            paths[itemp].mutate();
            paths[itemp].measure();

            if(itemp != ntemps-1) paths[itemp].temperature_swap_intra(paths[itemp+1]);
        }

        paths[0].temperature_swap_inter(true);          // start = true
        paths[ntemps-1].temperature_swap_inter(false);  // start = false
    }


    
    if(rank == 0) paths[0].print("best_path_angular_PT.out");

    double* l_process = new double[ntemps];
    for(int itemp = 0; itemp < ntemps; itemp++){
        l_process[itemp] = paths[itemp].get_length(1);
    }

    double* all_lengths = new double[size*ntemps];

    MPI_Gather(l_process,ntemps,MPI_REAL8,all_lengths,ntemps,MPI_REAL8,0,MPI_COMM_WORLD);
    if(rank == 0){
        ofstream Lengths("lengths_angular_PT.out");
        for(int i = 0; i < size*ntemps; i++) Lengths << temperatures[i] << '\t' << all_lengths[i] << endl;
        Lengths.close();

        cout << endl << "Length of the shortest path = " << paths[0].get_length(1) << endl << endl;
    }

    MPI_Finalize(); 
    return 0;
}
