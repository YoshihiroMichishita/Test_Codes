#include "C3v_HSL.hpp"

using namespace std;

void HSL::EigenV(parm parm_, int S){
    for (int i = 0; i < S; i++){
        //Complex H[M][M];
        double E0[M];
        Complex H[M*M];
        H_mom(parm_,KL[i],H);
        Diag_H<M>(H,E0);
        for (int j = 0; j < M; j++){
            E[i][j] = E0[j];
        }  
    }
}
/*
void HSL::Spectral(parm parm_, double re[][Mf], double im[][Mf], int S){
    for (int i = 0; i < S/5; i++){
        for(int WW=0; WW<parm_.W_SIZE; WW += 10){
            Complex G[M*M];
            GreenR_mom(parm_,KL[i],im[WW],re[WW],G);
            A[WW][i] = -imag(Trace_NH<M>(G));
        }
    }
}*/

//High Symmetric line
HSL::HSL(parm parm_){
    
    SIZE=0;
    cout << "HSL_.start" <<endl;
    double dk = 2 * pi / parm_.K_SIZE;
    double FF = 2.0/sqrt(3.0);
    //K => M => K'
    for (double KK = 0; KK < FF*2*pi/sqrt(3.0); KK += dk){
        double k0[D];
        k0[0] = FF*2*pi/sqrt(3.0) - KK/2; k0[1] = sqrt(3.0)*KK/2;; 
        double q = KK; 
        for (int i = 0; i < D; i++){
            KL[SIZE][i] = k0[i];
        }
        QL[SIZE] =q;
        
        SIZE++;
    }

    //K' => M
    for (double KM = 0; KM < FF*pi/sqrt(3.0); KM += dk){
        double k0[D];
        k0[0] = FF*pi/sqrt(3.0) - KM; k0[1] = FF*pi;
        double q = FF*2*pi/sqrt(3.0) + KM; 
        for (int i = 0; i < D; i++){
            KL[SIZE][i] = k0[i];
        }
        // length of symmetric line
        QL[SIZE] =q;
        
        SIZE++;
    }

    //M => G
    for (double MG = 0; MG < FF*pi; MG += dk){
        double k0[D];
        k0[0] = 0; k0[1] = FF*pi - MG;
        double q = FF*sqrt(3.0)*pi + MG; 
        for (int i = 0; i < D; i++){
            KL[SIZE][i] = k0[i];
        }
        // length of symmetric line
        QL[SIZE] =q;
        
        SIZE++;
    }

    //G => K
    for (double GK = 0; GK < FF*2*pi/sqrt(3.0); GK += dk){
        double k0[D];
        k0[0] = GK; k0[1] = 0; 
        double q = FF*(1.0 + sqrt(3.0)) * pi + GK; 
        for (int i = 0; i < D; i++){
            KL[SIZE][i] = k0[i];
        }
        QL[SIZE] =q;
        
        SIZE++;
    }
    /*
    for (double GK = 0; GK < pi; GK += dk){
        double k0[D];
        k0[0] = sqrt(3.0) * GK /2 ; k0[1] = GK/2;
        double q = (3.0 + sqrt(3.0)) * pi /2 + GK; 
        for (int i = 0; i < D; i++){
            KL[SIZE][i] = k0[i];
        }
        QL[SIZE] =q;
        
        SIZE++;
    }*/

    HSL::EigenV(parm_,SIZE);
     
};

HSL::~HSL(){
}
