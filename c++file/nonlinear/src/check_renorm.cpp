#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]){

    parm parm_(argv);
    parm_.Parm_List();

    ofstream out("Dirac_NHDisp.dat");

    
    double kx,ky;
    double w,im[2]={0,0},re[2]={0,0};
    cout << "kx ky w" << endl;
    cin >> kx >> ky >> w;
    double k[2] = {kx,ky};
    Ham Ham_(parm_,k,0);

    cout << "Ham" << endl;
    cout << Ham_.H_k[0] << "\t" << Ham_.H_k[1] << "\t" << Ham_.H_k[2] << "\t" << Ham_.H_k[3] << endl;
    cout << Ham_.VL_b[0] << "\t" << Ham_.VL_b[1] << "\t" << Ham_.VL_b[2] << "\t" << Ham_.VL_b[3] << endl;
    cout << Ham_.NH_fac[0] << "\t" << Ham_.NH_fac[1] << endl;
    cout << Ham_.VX_LL[0] << "\t" << Ham_.VX_LL[1] << "\t" << Ham_.VX_LL[2] << "\t" << Ham_.VX_LL[3] << endl;
    cout << Ham_.VX_RR[0] << "\t" << Ham_.VX_RR[1] << "\t" << Ham_.VX_RR[2] << "\t" << Ham_.VX_RR[3] << endl;

    Complex D1[M*M];
    Prod3<M>(Ham_.VL_b,Ham_.H_k,Ham_.VR_k,D1);
    cout << D1[0] << "\t" << D1[1] << "\t" << D1[2] << "\t" << D1[3] << endl;

    /*
    Green Green_(parm_,w,im,re,Ham_.H_k);
    cout << "Green" << endl;
    cout << Green_.GR[0] << "\t" << Green_.GR[1] << "\t" << Green_.GR[2] << "\t" << Green_.GR[3] << endl;
    Complex Gc[M*M],VRR[M*M],VLL[M*M],EE[M];
    Complex Gc2[M*M],VRR2[M*M],VLL2[M*M],EE2[M];
    for (int i = 0; i < M*M; i++){
        Gc[i] = Green_.GR[i];
        Gc2[i] = Green_.GA[i];
    }
    Diag_NH<M>(Gc,VLL,VRR,EE);
    Complex T0[M*M],T1[M*M],T2[M*M];
    Prod2<M>(VLL,VRR,T2);
    Prod3<M>(Ham_.VR_b,Green_.GA,Ham_.VL_k,T0);
    Prod3<M>(Ham_.VL_b,Green_.GR,Ham_.VR_k,T1);
    cout << "Test ";
    for (int i = 0; i < M*M; i++){
        cout << "\t" << T0[i];
        cout << "\t" << T1[i];
        //cout << "\t" << T2[i];
    }
    cout << endl;
    cout << endl;
    
    Diag_NH<M>(Gc2,VLL2,VRR2,EE2);
    
    double BI = 0;
    Complex GRR[M],GAA[M];
    for (int i = 0; i < M; i++){
        GRR[i] = 1.0/(w-Ham_.E_NH[i]+I*parm_.delta);
        GAA[i] = conj(GRR[i]);
    }
    cout << "compare (BI,Green) " << endl;
    cout << "Eigenvalue " << GRR[0] << "\t" << GRR[1] << "\t" << EE[0] << "\t" << EE[1] << endl;
    cout << "EigenvalueA " << GAA[0] << "\t" << GAA[1] << "\t" << EE2[0] << "\t" << EE2[1] << endl;
    cout << "VL ";
    for (int i = 0; i < M*M; i++){
        cout << "\t" << Ham_.VL_b[i];
        cout << "\t" << VLL[i];
    }
    cout << endl;
    cout << endl;
    cout << "VR ";
    for (int i = 0; i < M*M; i++){
        cout << "\t" << Ham_.VR_k[i];
        cout << "\t" << VRR[i];
    }
    cout << endl;
    cout << endl;
    cout << "check renorm " << endl;
    Complex U_LR[M*M],U_RL[M*M];
    Prod2<M>(Ham_.VL_k,Ham_.VR_b,U_LR);
    Prod2<M>(Ham_.VR_k,Ham_.VL_b,U_RL);
    for (int i = 0; i < M*M; i++){
        cout << "\t" << U_LR[i];
        cout << "\t" << U_RL[i];
    }
    cout << endl;


    for (int i = 0; i < M; i++){
        BI += real(Ham_.VX_RR[i*M+i] * GRR[i] * Ham_.VX_LL[i*M+i] * GAA[i]);
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                BI += real(Ham_.VX_RR[j*M+i] * GRR[i] * Ham_.VX_LL[i*M+j] * GAA[j]);
            }
            
        }
    }
    Complex GR2[M*M],GA2[M*M];
    for (int i = 0; i < M; i++){
        GR2[i*(M+1)] = GRR[i];
        GA2[i*(M+1)] = GAA[i];
    }
    Complex VXG2[M*M],VXGA2[M*M],MXX2[M*M];
    Prod2<M>(Ham_.VX_RR,GR2,VXG2);
    Prod2<M>(Ham_.VX_LL,GA2,VXGA2);
    Prod2<M>(VXG2,VXGA2,MXX2);
    double BI2 = Trace_H<M>(MXX2);
    cout << "==================" << endl;
    cout << "BI " << BI << "\t" << BI2 << endl;
    Complex VXG[M*M],VXGA[M*M],MXX[M*M];
    Prod2<M>(Ham_.VX,Green_.GR,VXG);
    Prod2<M>(Ham_.VX,Green_.GA,VXGA);
    Prod2<M>(VXG,VXGA,MXX);
    double GG = Trace_H<M>(MXX);
    cout << "GG " << GG << "\t" << MXX[0] << "\t" << MXX[3] << endl;*/

    /*for (int i = 0; i < parm_.K_SIZE; i++){
        double kx = 0.25 * i * pi /parm_.K_SIZE + 0.25*pi; 
        cout << "==================" << endl;
        cout << "kx " << kx << endl;
        double k[2] = {kx,0};
        Ham Ham_(parm_,k,2);
        out << kx << "\t" << real(Ham_.E_NH[0]) << "\t" << real(Ham_.E_NH[1]) << "\t" << imag(Ham_.E_NH[0]) << "\t" << imag(Ham_.E_NH[1]) << endl;
        Complex U0[M*M];
        Prod2<M>(Ham_.VL_b,Ham_.VR_k,U0);
        cout << "good basis ?" << endl;
        cout << U0[0] << "\t" << U0[1] << "\t" << U0[2] << "\t" << U0[3] << endl;
        cout << endl;
        cout << "Eigen Value" << endl;
        cout << Ham_.E_NH[0] << "\t" << Ham_.E_NH[1] << endl;

        Complex E0[M*M],RR[M*M],LL[M*M],E1[M];
        Prod3<M>(Ham_.VL_b,Ham_.H_k,Ham_.VR_k,E0);
        
        cout << E0[0] << "\t" << E0[1] << "\t" << E0[2] << "\t" << E0[3] << endl;
        cout << endl;
        cout << "NH factor" << endl;
        cout << Ham_.NH_fac[0] << "\t" << Ham_.NH_fac[1] << endl;
    }*/

    return 0;
}