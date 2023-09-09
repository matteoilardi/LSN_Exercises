#include "travelling_salesman.h"

using namespace std;

Map::Map(string fileName, Random* pr){
    m_pr = pr;

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

Map::Map(int ncities, int type, Random* pr){
    m_ncities = ncities;
    m_pr = pr;
        
    // circle
    if(type == 0){
        double theta;
        for(int i = 0; i<ncities; i++){
            theta = m_pr->Rannyu(0., 2*M_PI);
            m_x.push_back(cos(theta));
            m_y.push_back(sin(theta));
        }

        typeName = "circle";
    }
    // square
    else{
        for(int i = 0; i<ncities; i++){
            m_x.push_back(m_pr->Rannyu());
            m_y.push_back(m_pr->Rannyu());
        }

        typeName = "square";
    }

    dist_matrix.resize(ncities);
    for(int i = 0; i<ncities; i++){
        dist_matrix[i].resize(ncities);
        for(int j = 0; j<ncities; j++){
            dist_matrix[i][j] = sqrt(pow(m_x[i]-m_x[j], 2) + pow(m_y[i]-m_y[j], 2));
        }
    }
}


Path::Path(Map* pmap, Random* pr){
    m_pmap = pmap;
    m_pr = pr;
    m_ncities = m_pmap->get_ncities();

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
    Out.open(fileName+".out");

    for(auto icity: m_path) Out << m_pmap->get_x(icity) << '\t' << m_pmap->get_y(icity) << endl;

    Out.close();
}

string Path::get_map_type(){
    return m_pmap->get_type();
}

void Path::permutation(){
    int icity = int(m_pr->Rannyu(1., m_ncities));
    int jcity = pbc(icity + loaded_dice());

    swap(m_path[icity], m_path[jcity]);
}

void Path::shift(){
    int icity = int(m_pr->Rannyu(1., m_ncities));
    int n = loaded_dice(); // number of shifted cities
    int m = loaded_dice(); // shift

    for(int im = 0; im < m; im++){
        for(int in = n-1; in >= 0; in--){
            swap(m_path[pbc(icity+in+im)], m_path[pbc(icity+in+im + 1)]);
        }
    }
}

void Path::block_swap(){
    int icity = int(m_pr->Rannyu(1., m_ncities));
    int m = loaded_dice()/2; // block length

    for(int im = 0; im < m; im++){
        swap(m_path[pbc(icity+im)], m_path[pbc(icity+im + m)]);
    }
}

void Path::inversion(){
    int icity = int(m_pr->Rannyu(1., m_ncities));
    int m = loaded_dice();

    for(int im = 0; im < m/2; im++){
        swap(m_path[pbc(icity + im)], m_path[pbc(icity+m-1 - im)]);
    }
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
    double r = m_pr->Rannyu();
    int m = 1 + int((m_ncities-1) * pow(r, p));
    return m;
}



Population::Population(Map& map, int npaths, Random* pr){
    m_npaths = npaths;
    m_ncities = map.get_ncities();
    m_pr = pr;

    for(int ipath = 0; ipath < m_npaths; ipath++){
        Path p(&map, m_pr);
        m_pop.push_back(p);
    }
}

void Population::mutate(){
    for(int ipath = 0; ipath < m_npaths; ipath++){
        for(int icity = 0; icity < m_ncities; icity++){
            if(m_pr->Rannyu()<P2) m_pop[ipath].permutation();
        }

        if(m_pr->Rannyu()<P3) m_pop[ipath].shift();

        if(m_pr->Rannyu()<P4) m_pop[ipath].block_swap();

        if(m_pr->Rannyu()<P5) m_pop[ipath].inversion();
    }
}

void Population::newgen(){
    double r, s;
    int jpath, kpath;
    for(int ipath = 0; ipath < m_npaths/2; ipath++){
        r = m_pr->Rannyu();
        s = m_pr->Rannyu();
        jpath = int(m_npaths * pow(r, P6));
        kpath = int(m_npaths * pow(s, P6));

        crossover(jpath, kpath);
    }

    sort_paths();
    for(int ipath = 0; ipath < m_npaths; ipath++) m_pop.pop_back();
    // m_pop.resize(m_npaths);
}

void Population::crossover(int ipath, int jpath){
    int cut = int(m_pr->Rannyu(1., (double)m_ncities-2.));

    Path newborn1(m_pop[ipath]);
    Path newborn2(m_pop[jpath]);
    newborn1.get_path().resize(cut);
    newborn2.get_path().resize(cut);

    for(int kcity = 0; kcity < m_ncities; kcity++){
        if(! newborn1.is_in(m_pop[jpath].get_city(kcity))) newborn1.get_path().push_back(m_pop[jpath].get_city(kcity));
        if(! newborn2.is_in(m_pop[ipath].get_city(kcity))) newborn2.get_path().push_back(m_pop[ipath].get_city(kcity));
    }

    m_pop.push_back(newborn1);
    m_pop.push_back(newborn2);
}

void Population::sort_paths(){
    for(int ipath = 0; ipath < m_npaths; ipath++) m_pop[ipath].measure();
    sort(m_pop.begin(), m_pop.end());
}

double Population::half_av_l1(){
    double average = 0;
    for(int ipath = 0; ipath < m_npaths/2; ipath++){
        average += m_pop[ipath].get_length(1);
    }
    average = average/(m_npaths/2.);
    return average;
}

void Population::print_half_av_l1(int igen){
    ofstream Av;
    Av.open("half_av_"+m_pop[0].get_map_type()+".out", ios::app);

    Av << igen << '\t' << half_av_l1() << endl;
    Av.close();
}

Path Population::get_best_path_l1(){
    sort_paths();
    return m_pop[0];
}
