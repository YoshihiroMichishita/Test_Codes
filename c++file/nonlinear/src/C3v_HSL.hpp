#pragma once

#include "Ham_TMD.hpp"

using namespace std;

//This is the class for checking the dispersion on the high-symmetric line 

class HSL
{
    public:
    double KL[3000][D];
    double QL[3000];
    double E[3000][M];
    //double A[500][1000];
    int SIZE;

    HSL(parm parm_);
    void EigenV(parm parm_, int S);
    //void Spectral(parm parm_, int S);

    ~HSL();
};