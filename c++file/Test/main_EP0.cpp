#include "matrix_op_mypc.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
//#include <omp.h>

using namespace std;
typedef std::complex<double> Complex;

double t;
double Gamma;
double eta;

void Ham0(double k[2], Complex E[2]){
    Complex H[4], VL[4], VR[4];
    H[0] = t*k[0]+I*eta; H[1] = t*k[1];
    H[2] = t*k[1]; H[0] = t*k[0]+I*eta;
    Complex E0[2];
    Diag_NH<2>(H,VL,VR,E0);
    if(real(E0[0])>real(E0[1])){
        E[0] = E0[1]; E[1] = E0[0];
    }
    else{
        E[0] = E0[0]; E[1] = E0[1];
    }
}

void Ham1(double k[2], Complex E[2]){
    Complex H[4], VL[4], VR[4];
    H[0] = t*k[0]+I*(eta+Gamma); H[1] = t*k[1];
    H[2] = t*k[1]; H[0] = -t*k[0]+I*(eta-Gamma);
    Complex E0[2];
    Diag_NH<2>(H,VL,VR,E0);
    if(real(E0[0])>real(E0[1])){
        E[0] = E0[1]; E[1] = E0[0];
    }
    else{
        E[0] = E0[0]; E[1] = E0[1];
    }
}

void Ham2(double k[2], Complex E[2]){
    Complex H[4], VL[4], VR[4];
    H[0] = I*(eta+Gamma); H[1] = t*(k[1]-I*k[0]);
    H[2] = t*(k[1]+I*k[0]); H[0] = I*(eta-Gamma);
    Complex E0[2];
    Diag_NH<2>(H,VL,VR,E0);
    if(real(E0[0])>real(E0[1])){
        E[0] = E0[1]; E[1] = E0[0];
    }
    else{
        E[0] = E0[0]; E[1] = E0[1];
    }
}

int main(int argc, char* argv[]){

    int K_SIZE = atoi(argv[1]);
    t = atof(argv[2]);
    Gamma = atof(argv[3]);
    eta = atof(argv[4]);
    double K_MAX = atof(argv[5]);

    ofstream out("Dirac.dat");
    ofstream out1("normal_EP.dat");
    ofstream out2("PT_EP.dat");

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    for (int i = 0; i < K_SIZE; i++){
        for (int j = 0; j < K_SIZE; j++){
            double k[2] = {(i-K_SIZE/2)*K_MAX/K_SIZE, (j-K_SIZE/2)*K_MAX/K_SIZE};
            Complex E0[2],E1[2],E2[2];
            Ham0(k,E0);
            Ham1(k,E1);
            Ham2(k,E2);
            out << k[0] << "\t" << k[1] << "\t" << real(E0[0]) << "\t" << real(E0[1]) << endl;
            out1 << k[0] << "\t" << k[1] << "\t" << real(E1[0]) << "\t" << real(E1[1]) << endl;
            out2 << k[0] << "\t" << k[1] << "\t" << real(E2[0]) << "\t" << real(E2[1]) << endl;

        }
        out << endl;
        out1 << endl;
        out2<< endl;
    }

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);
    
    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    
    return 0;
}
