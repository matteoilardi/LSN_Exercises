#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;


class RandomWalk{
    public:
    RandomWalk(double step, Random* prand){
        m_step = step;
        m_start = {0., 0., 0.};
        m_r = m_start;
        m_prand = prand;
    };
    RandomWalk(double step, double x_start, double y_start, double z_start, Random* prand){
        m_step = step;
        m_start = {x_start, y_start, z_start};
        m_r = m_start;
        m_prand = prand;
    };
    virtual ~RandomWalk(){;};

    virtual void step() =  0;
    double r2(){
        double r2 = 0;
        for(auto x: m_r) r2 += x*x;
        return r2;
    };


    protected:

    double m_step;
    vector<double> m_start;
    vector<double> m_r;

    Random* m_prand;
};

class RW_lattice: public RandomWalk{
    public:
    RW_lattice(double step, Random* prand): RandomWalk(step, prand){;};
    RW_lattice(double step, double x_start, double y_start, double z_start, Random* prand): RandomWalk(step, x_start, y_start, z_start, prand){;};

    void step() override{
        int i = m_prand->randInt(0, 2);
        bool b = m_prand->randBool();

        b? m_r[i]+=m_step : m_r[i]-=m_step;
    };
};

class RW_continuum: public RandomWalk{
    public:
    RW_continuum(double step, Random* prand): RandomWalk(step, prand){;};
    RW_continuum(double step, double x_start, double y_start, double z_start, Random* prand): RandomWalk(step, x_start, y_start, z_start, prand){;};

    void step() override{
        double theta = m_prand->rand_theta();
        double phi = m_prand->rand_phi();

        m_r[0] += m_step*sin(theta)*cos(phi);
        m_r[1] += m_step*sin(theta)*sin(phi);
        m_r[2] += m_step*cos(theta);
    };
};
