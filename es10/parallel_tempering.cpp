#include "parallel_tempering.h"

using namespace std;


Map::Map(string fileName){

    ifstream input(fileName);
    double longit, latit;
    // const double R = 6378;
    m_ncities = 0;
    
    input >> longit >> latit;
    while(!input.eof()){
        // m_x.push_back(longit*M_PI/180. * R * sin((90.-latit)*M_PI/180.));
        // m_y.push_back(latit*M_PI/180. * R);
        m_x.push_back(longit);
        m_y.push_back(latit);
        m_ncities++;

        input >> longit >> latit;
    }
    input.close();


    dist_matrix.resize(m_ncities);
    for(int i = 0; i < m_ncities; i++){
        dist_matrix[i].resize(m_ncities);
        for(int j = 0; j < m_ncities; j++){
            dist_matrix[i][j] = sqrt(pow(m_x[i]-m_x[j], 2) + pow(m_y[i]-m_y[j], 2));
        }
    }
}


Path::Path(Map* pmap, double T, int rank, int size){
    int seed[4];
    int p1, p2, ptemp;
    ifstream Primes("Primes");
    if(Primes.is_open()){
        for(int irank = 0; irank < rank; irank++){
            Primes >> ptemp >> ptemp;
        }
        Primes >> p1 >> p2;
    }else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if(input.is_open()){
        while (!input.eof()){
            input >> property;
            if(property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                m_rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }else cerr << "PROBLEM: Unable to open seed.in" << endl;

    m_pmap = pmap;
    m_ncities = m_pmap->get_ncities();

    m_T = T;
    m_rank = rank;
    m_size = size;

    for(int i = 0; i < m_ncities; i++) m_path.push_back(i);
    for(int i = 0; i < m_ncities*100; i++) permutation();
}

bool Path::check(){
    vector<int> count(m_ncities, 0);
    bool ok = true;
        
    for(int i = 0; i < m_ncities; i++){
        count[m_path[i]]++;
    }
    for(int i = 0; i < m_ncities && ok == true; i++){
        if(count[i] != 1) ok = false;
    }

    return ok;
}
    
void Path::measure(){
    double length1 = 0.;
    double length2 = 0.;
    for(int icity = 0; icity < m_ncities-1; icity++){
        length1 += m_pmap->distance1(m_path[icity+1], m_path[icity]);
        length2 += m_pmap->distance2(m_path[icity+1], m_path[icity]);
    }
    length1 += m_pmap->distance1(m_path[0], m_path[m_ncities-1]);
    length2 += m_pmap->distance2(m_path[0], m_path[m_ncities-1]);

    m_length1 = length1;
    m_length2 = length2;
}

void Path::print(string fileName){
    ofstream Out;
    Out.open(fileName);

    for(auto icity: m_path) Out << m_pmap->get_x(icity) << '\t' << m_pmap->get_y(icity) << endl;

    Out.close();
}


void Path::permutation(){
    int icity = int(m_rnd.Rannyu(1., m_ncities));
    int jcity = pbc(icity + loaded_dice());

    swap(m_path[icity], m_path[jcity]);
}

void Path::shift(){
    int icity = int(m_rnd.Rannyu(1., m_ncities));
    int n = loaded_dice(); // number of shifted cities
    int m = loaded_dice(); // shift

    for(int im = 0; im < m; im++){
        for(int in = n-1; in >= 0; in--){
            swap(m_path[pbc(icity+in+im)], m_path[pbc(icity+in+im + 1)]);
        }
    }
}

void Path::block_swap(){
    int icity = int(m_rnd.Rannyu(1., m_ncities));
    int m = loaded_dice()/2; // block length

    for(int im = 0; im < m; im++){
        swap(m_path[pbc(icity+im)], m_path[pbc(icity+im + m)]);
    }
}

void Path::inversion(){
    int icity = int(m_rnd.Rannyu(1., m_ncities));
    int m = loaded_dice();

    for(int im = 0; im < m/2; im++){
        swap(m_path[pbc(icity + im)], m_path[pbc(icity+m-1 - im)]);
    }
}


void Path::temperature_swap_intra(Path& path_up){
    double length_up = path_up.get_length(1);
    double T_up = path_up.get_temperature();

    if(m_rnd.Rannyu() < p_swap_T_Boltzmann(length_up, T_up)){
        int* new_path_down = path_up.to_array();
        int* new_path_up = to_array();

        set_path(new_path_down);
        path_up.set_path(new_path_up);

        measure();
        path_up.measure();
    }
}

void Path::temperature_swap_inter(bool start){
    int* trial_path = new int[m_ncities];
    double trial_length;
    double trial_T;
    double T = get_temperature();
    int* path = to_array();
    double length = get_length(1);
    int accept;

    int itag_walker = 0;
    MPI_Status stat1, stat2, stat3, stat4, stat5;

    measure();

    if(m_rank != 0 and start == true){
        MPI_Send(path,m_ncities,MPI_INTEGER,m_rank-1,itag_walker+m_rank/2,MPI_COMM_WORLD);
        itag_walker += m_size/2;
        MPI_Send(&T,1,MPI_REAL8,m_rank-1,itag_walker+m_rank/2,MPI_COMM_WORLD);
        itag_walker += m_size/2;
        MPI_Send(&length,1,MPI_REAL8,m_rank-1,itag_walker+m_rank/2,MPI_COMM_WORLD);
        itag_walker += m_size/2;

        MPI_Recv(&accept,1,MPI_INTEGER,m_rank-1,itag_walker+m_rank/2, MPI_COMM_WORLD, &stat4);
        itag_walker += m_size/2; 

        if(accept){
            MPI_Recv(trial_path,m_ncities,MPI_INTEGER,m_rank-1,itag_walker+m_rank/2, MPI_COMM_WORLD, &stat5);
            itag_walker += m_size/2;
            set_path(trial_path);
        }

    }else if(m_rank != m_size-1 and start == false){
        MPI_Recv(trial_path,m_ncities,MPI_INTEGER,m_rank+1,itag_walker+(m_rank+1)/2, MPI_COMM_WORLD, &stat1);
        itag_walker += m_size/2;
        MPI_Recv(&trial_T,1,MPI_REAL8,m_rank+1,itag_walker+(m_rank+1)/2, MPI_COMM_WORLD, &stat2);
        itag_walker += m_size/2;
        MPI_Recv(&trial_length,1,MPI_REAL8,m_rank+1,itag_walker+(m_rank+1)/2, MPI_COMM_WORLD, &stat3);
        itag_walker += m_size/2;

        if(m_rnd.Rannyu() < p_swap_T_Boltzmann(trial_length, trial_T)) accept = 1;
        else accept = 0;
        MPI_Send(&accept,1,MPI_INTEGER,m_rank+1,itag_walker+(m_rank+1)/2,MPI_COMM_WORLD);
        itag_walker += m_size/2;

        if(accept){
            MPI_Send(trial_path,m_ncities,MPI_INTEGER,m_rank+1,itag_walker+(m_rank+1)/2,MPI_COMM_WORLD);
            itag_walker += m_size/2;
            set_path(trial_path);
        }
    }

    measure();
    delete[] trial_path;
}


void Path::mutate(){
    Path trial_path(*this);
    for(int icity = 0; icity < m_ncities; icity++){
        trial_path = *this;
        trial_path.permutation();
        if(m_rnd.Rannyu() < p_Boltzmann(trial_path)) *this = trial_path;
    }

    trial_path = *this;
    trial_path.shift();
    if(m_rnd.Rannyu() < p_Boltzmann(trial_path)) *this = trial_path;

    trial_path = *this;
    trial_path.block_swap();
    if(m_rnd.Rannyu() < p_Boltzmann(trial_path)) *this = trial_path;

    trial_path = *this;
    trial_path.inversion();
    if(m_rnd.Rannyu() < p_Boltzmann(trial_path)) *this = trial_path;
}


double Path::p_Boltzmann(Path trial_path){
    measure();
    trial_path.measure();
    double delta_length = trial_path.get_length(1) - get_length(1);
    return exp(-1./m_T*delta_length);
}


double Path::p_swap_T_Boltzmann(double trial_length, double trial_T){
    double delta_beta = 1./get_temperature() - 1./trial_T;
    double delta_length = trial_length - get_length(1);
    return exp(-1.*delta_beta*delta_length);
}



int Path::pbc(int icity){
    if(icity == 0) return 1;
    while(icity >= m_ncities){
        icity = icity - (m_ncities-1);
    }
    return icity;
}

bool Path::is_in(int number){
    for(auto city: m_path){
        if(city == number) return true;
    }
    return false;
}

int Path::loaded_dice(double p){
    double r = m_rnd.Rannyu();
    int m = 1 + int((m_ncities-1) * pow(r, p));
    return m;
}


vector<double> readTemperatures(string fileName){
    ifstream input(fileName);
    double T;
    vector<double> temperatures;

    input >> T;
    while(!input.eof()){
        temperatures.push_back(T);
        input >> T;
    }

    input.close();
    return temperatures;
}
