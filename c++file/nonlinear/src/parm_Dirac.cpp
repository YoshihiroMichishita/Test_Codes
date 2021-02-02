#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "const.hpp"
#include "parm_Dirac.hpp"



using namespace std;

parm::parm(char* argv[]){

    K_SIZE = atoi(argv[1]);
    t_i = atof(argv[2]);
    t_e = atof(argv[3]);
    delta = atof(argv[4]);
    alpha = atof(argv[5]);
    mu = atof(argv[6]);
    Ma = atof(argv[7]);
    T = atof(argv[8]);
    NH = atof(argv[9]);
    TB = atof(argv[10]);
    W_SIZE = atoi(argv[11]);
    W_MAX = atof(argv[12]);
    W = atof(argv[13]);
    K_MAX = atof(argv[14]);
  
}


parm::~parm(){  
}

void parm::Parm_List(){
    cout<<"********************"<<endl;
    cout<<"K_SIZE "<<K_SIZE<<endl;
    cout<<"K_MAX "<<K_MAX<<endl;
    cout<<"t_x,y "<<t_i<<endl;
    cout<<"gamma "<<t_e<<endl;
    cout<<"delta "<<delta<<endl;
    cout<<"alpha "<<alpha<<endl;
    cout<<"mu "<<mu<<endl;
    cout<<"Mass "<<Ma<<endl;
    cout<<"T "<<T<<endl;
    cout<<"NH "<<NH<<endl;
    cout<<"TB "<<TB<<endl;
    cout<<"W "<<W<<endl;
    cout<<"W_SIZE "<<W_SIZE<<endl;
    cout<<"W_MAX "<<W_MAX<<endl;
    cout<<"********************"<<endl;
}