#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include "const.hpp"
#include "parm_sDirac.hpp"



using namespace std;

parm::parm(char* argv[]){

    K_SIZE = atoi(argv[1]);
    vx = atof(argv[2]);
    vy = atof(argv[3]);
    delta = atof(argv[4]);
    alpha = atof(argv[5]);
    im_b = atof(argv[6]);
    re_b = atof(argv[7]);
    mu = atof(argv[8]);
    T = atof(argv[9]);
    W_MAX = atof(argv[10]);
    W_SIZE = atof(argv[11]);
    K_MAX = atof(argv[12]);
  
}


parm::~parm(){  
}

void parm::Parm_List(){
    cout<<"********************"<<endl;
    cout<<"K_SIZE "<<K_SIZE<<endl;
    cout<<"K_MAX "<<K_MAX<<endl;
    cout<<"vx "<<vx<<endl;
    cout<<"vy "<<vy<<endl;
    cout<<"delta "<<delta<<endl;
    cout<<"alpha "<<alpha<<endl;
    cout<<"beta(im,re)  (" <<im_b << ", " << re_b << ") " << endl;
    cout<<"mu "<<mu<<endl;
    cout<<"T "<<T<<endl;
    cout<<"W_SIZE "<<W_SIZE<<endl;
    cout<<"W_MAX "<<W_MAX<<endl;
    cout<<"********************"<<endl;
}