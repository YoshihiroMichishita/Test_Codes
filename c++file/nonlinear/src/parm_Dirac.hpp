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

    double t_i;
    double t_e;
    double delta;
    double alpha;
    double mu;
    double Ma;
    double W;
    double T;
    double W_MAX;
    int W_SIZE;
    double NH;
    double TB;

    parm(char *argv[]);
    ~parm();
    void Parm_List();
};