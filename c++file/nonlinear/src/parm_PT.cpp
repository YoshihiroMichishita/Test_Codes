#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "const.hpp"
#include "parm_PT.hpp"



using namespace std;

parm::parm(char* argv[]){

    K_SIZE = atoi(argv[1]);
    t_i = atof(argv[2]);
    t_e = atof(argv[3]);
    delta = atof(argv[4]);
    a_R = atof(argv[5]);
    a_D = atof(argv[6]);
    mu = atof(argv[7]);
    hx = atof(argv[8]);
    hy = atof(argv[9]);
    hz = atof(argv[10]);
    T = atof(argv[11]);
    NH = atof(argv[12]);
    TB = atof(argv[13]);
    W_SIZE = atoi(argv[14]);
    W_MAX = atof(argv[15]);
    W = atof(argv[16]);
  
}


parm::~parm(){  
}

void parm::Parm_List(){
    cout<<"********************"<<endl;
    cout<<"K_SIZE "<<K_SIZE<<endl;
    cout<<"t_i "<<t_i<<endl;
    cout<<"t_e "<<t_e<<endl;
    cout<<"delta "<<delta<<endl;
    cout<<"a_R, a_D "<<a_R << ", " << a_D <<endl;
    cout<<"mu "<<mu<<endl;
    cout<<"(hx,hy,hz) = (" << hx << ", " << hy << ", " << hz << ") " <<endl;
    cout<<"T "<<T<<endl;
    cout<<"NH "<<NH<<endl;
    cout<<"TB "<<TB<<endl;
    cout<<"W "<<W<<endl;
    cout<<"W_SIZE "<<W_SIZE<<endl;
    cout<<"W_MAX "<<W_MAX<<endl;
    cout<<"********************"<<endl;
}