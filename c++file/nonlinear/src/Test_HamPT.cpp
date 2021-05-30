#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;
typedef std::complex<double> Complex;

int main(int argc, char* argv[]){
    parm parm_(argv);
    parm_.Parm_List();
    cout << "kx,ky,w" << endl;
    double k[2];
    int sw=0;
    double w;
    cin >> k[1] >> k[2] >> w;

    double ff = FD(w,parm_.T);
    Ham Ham_(parm_,k,0);
    Ham_.Ham_list(0);
    double ddk = 4.0 * pi * pi /(parm_.K_SIZE*parm_.K_SIZE);
    double PVC2=0,PVC3=0,PVC=0;
    double G_xxy1=0,GBI_xxy1=0,G_xxy2=0,GBI_xxy2=0;
    double G_xxy3=0,GBI_xxy3=0,G_xxy4=0,GBI_xxy4=0;
    double re[2]={0,0},im[2]={0,0};


    Green Green_(parm_,w,im,re,Ham_.H_k);

    /*
    Complex VXXGRP[M*M],VYGRA[M*M],VXYGRM[M*M],VXGRA[M*M];
    Prod2<M>(Ham_.VXX,Green_.GRp,VXXGRP);
    Prod2<M>(Ham_.VY,Green_.GRmA,VYGRA);
    Prod2<M>(Ham_.VYX,Green_.GRm,VXYGRM);
    Prod2<M>(Ham_.VX,Green_.GRmA,VXGRA);

    Complex VXXGRA[M*M],VYGAM[M*M],VXYGRA[M*M],VXGAP[M*M];
    Prod2<M>(Ham_.VXX,Green_.GRmA,VXXGRA);
    Prod2<M>(Ham_.VY,Green_.GAm,VYGAM);
    Prod2<M>(Ham_.VYX,Green_.GRmA,VXYGRA);
    Prod2<M>(Ham_.VX,Green_.GAp,VXGAP);

    Complex PV3_1[M*M],PV3_2[M*M],PV3_3[M*M],PV3_4[M*M];
    Prod2<M>(VXXGRP,VYGRA,PV3_1);
    Prod2<M>(VXYGRM,VXGRA,PV3_2);

    Prod2<M>(VXXGRA,VYGAM,PV3_3);
    Prod2<M>(VXYGRA,VXGAP,PV3_4);
    G_xxy1 = Trace_H<M>(PV3_1);
    G_xxy2 = Trace_H<M>(PV3_2);
    G_xxy3 = Trace_H<M>(PV3_3);
    G_xxy4 = Trace_H<M>(PV3_4);
    //Opt_Circular_Green(parm_,Ham_,Green_,G_xxy1,PVC2,PVC3,PVC);

    double PV1 = 0;
    
    
    for (int i = 0; i < M; i++){
        Complex GRPi = 1.0/(w+parm_.W-Ham_.EN[i]+I*parm_.delta);
        Complex GRMi = 1.0/(w-parm_.W-Ham_.EN[i]+I*parm_.delta);
        Complex GRmAi = 1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta);
        for (int j = 0; j < M; j++){
            Complex GAPj = 1.0/(w+parm_.W-Ham_.EN[j]-I*parm_.delta);
            Complex GAMj = 1.0/(w-parm_.W-Ham_.EN[j]-I*parm_.delta);
            Complex GRmAj = I*imag(2.0/(w-Ham_.EN[j]+I*parm_.delta));
            GBI_xxy1 += real(Ham_.VXX_LR[j*M+i]*GRPi*Ham_.VY_LR[i*M+j]*GRmAj);
            GBI_xxy3 += real(Ham_.VXX_LR[j*M+i]*GRmAi*Ham_.VY_LR[i*M+j]*GAMj);
            GBI_xxy2 += real(Ham_.VYX_LR[j*M+i]*GRMi*Ham_.VX_LR[i*M+j]*GRmAj);
            GBI_xxy4 += real(Ham_.VYX_LR[j*M+i]*GRmAi*Ham_.VX_LR[i*M+j]*GAPj);
        }
    }*/
    Ham_.Ham_list(0);
    Opt_Circular_Green(parm_,Ham_,Green_,w,G_xxy2,G_xxy3,G_xxy4,G_xxy1);
    Ham_.Ham_list(0);
    Opt_Circular_Green_BI(parm_,Ham_,w,GBI_xxy2,GBI_xxy1);

    cout << G_xxy1 << " " << G_xxy2 << " " << GBI_xxy1 << " " << GBI_xxy2 << endl;
    cout << G_xxy3 << " " << G_xxy4 << " " << GBI_xxy3 << " " << GBI_xxy4 << endl;

    //double Green_xxy_re=0;
    //Opt_Gyration_BI3(parm_,Ham_,BI_xxy_re,BI_xxy_im);
    //double re[2]={0,0},im[2]={0,0};
    //Green Green_(parm_,0,im,re,Ham_.H_k);

    //Green_.G_List();
    //Ham_.Ham_list(sw);

    return 0;
    
}