#ifndef __MonteCarlo__
#define __MonteCarlo__

#include <cmath>
#include "random.h"

using namespace std;

class distribution{
    public:
    distribution(double min, double max);
    distribution(){;};

    virtual double Rand() = 0;
    virtual double Eval(double x) const = 0;
    
    protected:
    double m_min, m_max;
    Random rnd;
};

class uniform: public distribution{
    public:
    uniform(double min, double max): distribution(min, max){;};

    double Rand() override{
        return rnd.Rannyu(m_min, m_max);
    }
    double Eval(double x) const override{
        return 1./(m_max-m_min);
    }
};

class circle_quarter: public distribution{
    public:
    circle_quarter(double min, double max): distribution(min, max){
        m_r = m_max-m_min;
        m_k = 4./M_PI/(m_r*m_r);
    }

    double Rand() override{
        double x, y;
        do{
            x = rnd.Rannyu()*m_r;
            y = rnd.Rannyu()*m_r*m_k;
        }while(y > m_k*sqrt(m_r*m_r - x*x));

        return x + m_min;
    }
    double Eval(double x) const override{
        return m_k*sqrt(m_r*m_r - pow(x-m_min, 2));
    }

    private:
    double m_r, m_k;
};

class linear: public distribution{
    public:
    linear(double min, double max, double y_max): distribution(min, max){
        m_y_min = 2./(max-min) - y_max;
        m_y_max = y_max;

        m_m = (m_y_max-m_y_min)/(m_max-m_min);
        m_q = m_y_min - m_m*m_min;
    }

    double Rand() override{
       double y = rnd.Rannyu();
       double x = 1./m_m * (-1.*m_q + sqrt(m_q*m_q + pow(m_m*m_min, 2) + 2*m_m*m_q*m_min + 2*m_m*y));
       return x;
    }
    double Eval(double x) const override{
        double y = (x-m_min)*(m_y_max-m_y_min)/(m_max-m_min) + m_y_min;
        return y;
    }

    private:
    double m_y_min, m_y_max, m_m, m_q;
};



double MC(distribution* p, double (*f)(double), int n_random_block);

#endif