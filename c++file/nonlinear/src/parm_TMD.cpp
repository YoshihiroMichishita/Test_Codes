#include <iostream>
#include <stdlib.h>
#include<string>
#include<sstream>
#include<fstream>
#include "const.hpp"
#include "parm_TMD.hpp"



using namespace std;

parm::parm(char* argv[]){

    K_SIZE = atoi(argv[1]);
    t_i = atof(argv[2]);
    t_e = atof(argv[3]);
    a_u = atof(argv[4]);
    a_d = atof(argv[5]);
    delta = atof(argv[6]);
    alpha = atof(argv[7]);
    mu = atof(argv[8]);
    T = atof(argv[9]);
    Pr = atof(argv[10]);
    NH = atof(argv[11]);
    TB = atof(argv[12]);
    W_SIZE = atoi(argv[13]);
    W_MAX = atof(argv[14]);
    W = atof(argv[16]);
  
}


parm::~parm(){  
}

void parm::Parm_List(){
    cout<<"********************"<<endl;
    cout<<"K_SIZE "<<K_SIZE<<endl;
    //cout<<"K_MAX "<<K_MAX<<endl;
    cout<<"t_intra "<<t_i<<endl;
    cout<<"t_inter "<<t_e<<endl;
    cout<<"a_u "<<a_u<<endl;
    cout<<"a_d "<<a_d<<endl;
    cout<<"delta "<<delta<<endl;
    cout<<"alpha "<<alpha<<endl;
    cout<<"mu "<<mu<<endl;
    cout<<"T "<<T<<endl;
    cout<<"Pr "<<Pr<<endl;
    cout<<"NH "<<NH<<endl;
    cout<<"TB "<<TB<<endl;
    cout<<"W_SIZE "<<W_SIZE<<endl;
    cout<<"W_MAX "<<W_MAX<<endl;
    cout<<"********************"<<endl;
}