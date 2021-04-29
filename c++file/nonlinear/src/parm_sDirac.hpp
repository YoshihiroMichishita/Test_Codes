#pragma once
//#include <vector>
#include "const.hpp"
typedef std::complex<double> Complex;

using namespace std;

class parm
{
    public:
    int K_SIZE;
    double K_MAX;

    double vx;
    double vy;
    double alpha;
    double im_b;
    double re_b;
    double mu;
    double delta;
    double T;
    double W;
    double W_MAX;
    int W_SIZE;

    parm(char *argv[]);
    ~parm();
    void Parm_List();
};