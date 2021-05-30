#include "transport.hpp"

using namespace std;

//Fermi distribution function 
double FD(double w, double T){
    double f = 1.0/(1.0+exp(w/T));
    return f;
};

//w-derivative of teh Fermi distribution function
double dFD(double w, double T){
    double f;
    if(fabs(w/T)>10.0){
        f=0;
    }
    else{
        f = 1.0/(T*(1.0+exp(w/T))*(1.0+exp(-w/T)));
    }
    return f;
};


double ddFD(double w, double T){
    double f;
    if(fabs(w/T)>10.0){
        f=0;
    }
    else{
        f = 1.0/(T*T*(1.0+exp(w/T))*(1.0+exp(-w/T)))*(1.0/(1.0+exp(-w/T))-1.0/(1.0+exp(w/T)));
    }
    return f;
};

//spectral function A(k,w)
double Spectral(parm parm_,double w, double k[D], double re[Mf], double im[Mf]){
    Complex G[M*M],H[M*M];
    H_mom(parm_,k,H);
    GreenR_mom(parm_,w,im,re,H,G);
    double A = -imag(Trace_NH<M>(G))/pi;
    return A;
};

// Calculate linear transport and nonlnear Hall conductivity
void DC_transport(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2){
    Complex VXG[M*M],VYG[M*M],VYGG[M*M],VYGGG[M*M],VXGA[M*M];
    Prod2<M>(Ham_.VX,Green_.GR,VXG);
    Prod2<M>(Ham_.VY,Green_.GR,VYG);
    Prod2<M>(Ham_.VY,Green_.dGR,VYGG);
    Prod2<M>(Ham_.VX,Green_.GA,VXGA);

    Complex MXX[M*M],MYX[M*M],MYXX1[M*M],MYXX[M*M],MYXX2[M*M];
    Prod2<M>(VXG,VXGA,MXX);
    Prod2<M>(VYG,VXGA,MYX);
    Prod3<M>(VYGG,VXG,VXGA,MYXX);
    Prod3<M>(VYGG,Ham_.VXX,Green_.GA,MYXX2);

    XX += dw * Trace_H<M>(MXX) * dFD(w,T);
    YX += dw * Trace_NHi<M>(MYXX1) * dFD(w,T);
    YXXi += 2.0 * dw * Trace_NHi<M>(MYXX) * dFD(w,T);
    YXXr += 2.0 * dw * Trace_H<M>(MYXX) * dFD(w,T);
    YXXr2 += 2.0 * dw * Trace_H<M>(MYXX2) * dFD(w,T);
    YXXi2 += 2.0 * dw * Trace_NHi<M>(MYXX2) * dFD(w,T); 
};

// Calculate linear transport and non-reciprcal conductivity
void DC_transport_NRC(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& XXXr, double& XXXr2, double& XXXi, double& XXXi2){
    Complex VXG[M*M],VXGG[M*M],VYGGG[M*M],VXGA[M*M];
    Prod2<M>(Ham_.VX,Green_.GR,VXG);
    Prod2<M>(Ham_.VX,Green_.dGR,VXGG);
    Prod2<M>(Ham_.VX,Green_.GA,VXGA);

    Complex MXX[M*M],MXXX[M*M],MXXX2[M*M];
    Prod2<M>(VXG,VXGA,MXX);
    Prod3<M>(VXGG,VXG,VXGA,MXXX);
    Prod3<M>(VXGG,Ham_.VXX,Green_.GA,MXXX2);

    XX += dw * Trace_H<M>(MXX) * dFD(w,T);
    //YX += dw * Trace_NHi<M>(MYXX1) * dFD(w,T);
    XXXi += 2.0 * dw * Trace_NHi<M>(MXXX) * dFD(w,T);
    XXXr += 2.0 * dw * Trace_H<M>(MXXX) * dFD(w,T);
    XXXr2 += 2.0 * dw * Trace_H<M>(MXXX2) * dFD(w,T);
    XXXi2 += 2.0 * dw * Trace_NHi<M>(MXXX2) * dFD(w,T); 
};

// Calculate linear transport, nonlnear Hall conductivity, and non-reciprcal conductivity
void DC_NLH_NRC(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& YX, double& XXXr, double& XXXr2, double& XXXi, double& XXXi2, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2){
    Complex VXG[M*M],VYG[M*M],VYGG[M*M],VXGG[M*M],VXGA[M*M];
    Prod2<M>(Ham_.VX,Green_.GR,VXG);
    Prod2<M>(Ham_.VY,Green_.GR,VYG);
    Prod2<M>(Ham_.VX,Green_.dGR,VXGG);
    Prod2<M>(Ham_.VY,Green_.dGR,VYGG);
    Prod2<M>(Ham_.VX,Green_.GA,VXGA);

    Complex MXX[M*M],MYX[M*M],MYXX1[M*M],MYXX[M*M],MYXX2[M*M];
    Complex MXXX[M*M],MXXX2[M*M];
    Prod3<M>(VXGG,VXG,VXGA,MXXX);
    Prod3<M>(VXGG,Ham_.VXX,Green_.GA,MXXX2);
    Prod2<M>(VXG,VXGA,MXX);
    Prod2<M>(VYG,VXGA,MYX);
    Prod3<M>(VYGG,VXG,VXGA,MYXX);

    XX += dw * Trace_H<M>(MXX) * dFD(w,T);
    YX += dw * Trace_NHi<M>(MYXX1) * dFD(w,T);
    XXXi += 2.0 * dw * Trace_NHi<M>(MXXX) * dFD(w,T);
    XXXr += 2.0 * dw * Trace_H<M>(MXXX) * dFD(w,T);
    XXXr2 += 2.0 * dw * Trace_H<M>(MXXX2) * dFD(w,T);
    XXXi2 += 2.0 * dw * Trace_NHi<M>(MXXX2) * dFD(w,T); 
    YXXi += 2.0 * dw * Trace_NHi<M>(MYXX) * dFD(w,T);
    YXXr += 2.0 * dw * Trace_H<M>(MYXX) * dFD(w,T);
    YXXr2 += 2.0 * dw * Trace_H<M>(MYXX2) * dFD(w,T);
    YXXi2 += 2.0 * dw * Trace_NHi<M>(MYXX2) * dFD(w,T); 
};

// Calculate linear transport and nonlnear Hall conductivity with non-Hermitian band-index
void Linear_transport_NH(double dw, double w,double delta, double T,Ham Ham_, double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi,
 double& YXXi2, double& XX_, double& YX_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_){
    
    Complex GR[M],GA[M];
    for (int i = 0; i < M; i++){
        GR[i] = 1.0/(w - Ham_.E_NH[i]+I*delta);
        GA[i] = conj(GR[i]);
    }
    Complex YXX0=0,YXX0_=0,YXX1=0,YXX1_=0;

    for (int i = 0; i < M; i++){
        XX += dw * real(Ham_.VX_RR[i*(M+1)] * GR[i] * Ham_.VX_LL[i*(M+1)] * GA[i]) * dFD(w,T);
        if(fabs(Ham_.NH_fac[0])<1000.0){
            XX_ += dw * real( (Ham_.VR_b[M*i]*Ham_.VR_k[i] + Ham_.VR_b[M*i+1]*Ham_.VR_k[i+M]) * Ham_.VX_LR[i*(M+1)] * GR[i] * Ham_.VX_LR[i*(M+1)] * (Ham_.VL_b[M*i]*Ham_.VL_k[i] + Ham_.VL_b[M*i+1]*Ham_.VL_k[i+M]) * GA[i]) * dFD(w,T);

            for (int j = 0; j < M; j++){
                for (int l = 0; l < M; l++){
                    //all contribution
                    YX += -2.0 * dw * imag(Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j]) * dFD(w,T);
                    YX_ += -2.0 * dw * imag((Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]
                    * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LR[l*M+j] * (Ham_.VL_b[M*j]*Ham_.VL_k[j] + Ham_.VL_b[M*j+1]*Ham_.VL_k[M+j]) * GA[j]) * dFD(w,T);
                }
                YXX0 += -2.0*dw * Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+j] * GR[j] * Ham_.VX_LL[j*(M+1)] * GA[j] * dFD(w,T);
                YXX0_ += -2.0*dw * (Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]
                    * Ham_.VX_LR[M*i+j] * GR[j] * Ham_.VX_LR[j*(M+1)] * (Ham_.VL_b[M*j]*Ham_.VL_k[j] + Ham_.VL_b[M*j+1]*Ham_.VL_k[M+j]) * GA[j] * dFD(w,T);
                
                YXX1 += -2.0*dw * Ham_.VY_RR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LL[j+M*i] * GA[j] * dFD(w,T);
                YXX1_ += -2.0*dw * (Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]
                * Ham_.VXX_LR[M*i+j] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
            }
        }
        else{
            XX_ += dw * real(Ham_.VX_RR[i*(M+1)] * GR[i] * Ham_.VX_LL[i*(M+1)] * GA[i]) * dFD(w,T);

            for (int j = 0; j < M; j++){
                for (int l = 0; l < M; l++){
                    //all contribution
                    double QQ = -2.0 * dw * imag(Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j]) * dFD(w,T);
                    YX += QQ;
                    YX_ += QQ;
                }
                Complex PP1 = -2.0*dw * Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+j] * GR[j] * Ham_.VX_LL[j*(M+1)] * GA[j] * dFD(w,T);
                YXX0 += PP1;
                YXX0_ += PP1;
                Complex PP2 = -2.0*dw * Ham_.VY_RR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LL[j+M*i] * GA[j] * dFD(w,T);
                YXX1 += PP2;
                YXX1_ += PP2;
            }
        }
            
        
    }
    YXXi += imag(YXX0);
    YXXr += real(YXX0);
    YXXi_ += imag(YXX0_);
    YXXr_ += real(YXX0_);
    YXXi2 += imag(YXX1);
    YXXr2 += real(YXX1);
    YXXi2_ += imag(YXX1_);
    YXXr2_ += real(YXX1_);
};

// Calculate linear transport, nonlnear Hall conductivity, and non-reciprcal conductivity with non-Hermitian band-index
void Linear_transport_NHwithNRC(parm parm_,double dw, double w,double T,Ham Ham_, double& XX, double& XXX, double& XXX2, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2, double& XX_, double& XXX_, double& XXX2_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_){
    
    Complex GR[M],GA[M];
    for (int i = 0; i < M; i++){
        GR[i] = 1.0/(w - Ham_.E_NH[i]+I*parm_.delta);
        GA[i] = conj(GR[i]);
    }
    Complex YXX0=0,YXX0_=0,YXX1=0,YXX1_=0;
    Complex XXX0=0,XXX0_=0,XXX1=0,XXX1_=0;

    for (int i = 0; i < M; i++){
        
        XX_ += dw * real( (Ham_.VR_b[M*i]*Ham_.VR_k[i] + Ham_.VR_b[M*i+1]*Ham_.VR_k[i+M]) * Ham_.VX_LR[i*(M+1)] * GR[i] * Ham_.VX_LR[i*(M+1)] * (Ham_.VL_b[M*i]*Ham_.VL_k[i] + Ham_.VL_b[M*i+1]*Ham_.VL_k[i+M]) * GA[i]) * dFD(w,T);

        for (int j = 0; j < M; j++){
            XX += dw * real(Ham_.VX_RR[j*M+i] * GR[i] * Ham_.VX_LL[i*M+j] * GA[j]) * dFD(w,T);
            //XXX0 += dw * Ham_.VX_RR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[j+M*i] * GR[j] * Ham_.VX_LL[j*(M+1)] * GA[j] * dFD(w,T);
            //XXX0_ += dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[j+M*i] * GR[j] * Ham_.VX_LR[j*(M+1)] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
            
            XXX1 += -dw * Ham_.VX_RR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LL[j+M*i] * GA[j] * dFD(w,T);
            XXX1_ += -dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LR[j+M*i] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);

            for (int l = 0; l < M; l++){
                //all contribution
                XXX0 += -dw * Ham_.VX_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j] * dFD(w,T);
                XXX0_ += -dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LR[l*M+j] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
                if(j==l){

                }
                else{
                    YXX0 += -2.0*dw * Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j] * dFD(w,T);
                    YXX0_ += -2.0*dw * (Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]* Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LR[l*M+j] * (Ham_.VL_b[M*j]*Ham_.VL_k[j] + Ham_.VL_b[M*j+1]*Ham_.VL_k[M+j]) * GA[j] * dFD(w,T);
                }
                //YX += -2.0 * dw * imag(Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j]) * dFD(w,T);
                //YX_ += -2.0 * dw * imag((Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]
                // * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LR[l*M+j] * (Ham_.VL_b[M*j]*Ham_.VL_k[j] + Ham_.VL_b[M*j+1]*Ham_.VL_k[M+j]) * GA[j]) * dFD(w,T);
            }
            YXX0 += -2.0*dw * Ham_.VY_RR[j*M+i] * GR[i] * GR[i] * Ham_.VX_LR[M*i+j] * GR[j] * Ham_.VX_LL[j*(M+1)] * GA[j] * dFD(w,T);
            YXX0_ += -2.0*dw * (Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i]
                * Ham_.VX_LR[M*i+j] * GR[j] * Ham_.VX_LR[j*(M+1)] * (Ham_.VL_b[M*j]*Ham_.VL_k[j] + Ham_.VL_b[M*j+1]*Ham_.VL_k[M+j]) * GA[j] * dFD(w,T);
            
            YXX1 += -2.0*dw * Ham_.VY_RR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LL[j+M*i] * GA[j] * dFD(w,T);
            YXX1_ += -2.0*dw * (Ham_.VR_b[M*j]*Ham_.VR_k[j] + Ham_.VR_b[M*j+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[j*M+i] * GR[i] * GR[i] * Ham_.VXX_LR[M*i+j] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
        }
        
    }
    YXXi += imag(YXX0);
    YXXr += real(YXX0);
    YXXi_ += imag(YXX0_);
    YXXr_ += real(YXX0_);
    YXXi2 += imag(YXX1);
    YXXr2 += real(YXX1);
    YXXi2_ += imag(YXX1_);
    YXXr2_ += real(YXX1_);
    XXX += imag(XXX0);
    XXX2 += imag(XXX1);
    XXX_ += imag(XXX0_);
    XXX2_ += imag(XXX1_);
};

// Calculate linear transport and non-reciprcal conductivity with non-Hermitian band-index
void Linear_transport_NH_NRC(double dw, double w,double T,Ham Ham_, double& XX, double& XXX, double& XX_, double& XXX_){
    Complex GR[M],GA[M];
    for (int i = 0; i < M; i++){
        GR[i] = 1.0/(w - Ham_.E_NH[i]);
        GA[i] = conj(GR[i]);
    }
    Complex XXX0=0,XXX0_=0,XXX1=0,XXX1_=0,XXX2=0,XXX2_=0;
    for (int i = 0; i < M; i++){
        XX += dw * real(Ham_.VX_RR[i*(M+1)] * GR[i] * Ham_.VX_LL[i*(M+1)] * GA[i]) * dFD(w,T);
        XX_ += dw * real( (Ham_.VR_b[i*M]*Ham_.VR_k[i] + Ham_.VR_b[i*M+1]*Ham_.VR_k[i+M]) * Ham_.VX_LR[i*(M+1)] * GR[i] * Ham_.VX_LR[i*(M+1)] * (Ham_.VL_b[i*M]*Ham_.VL_k[i] + Ham_.VL_b[i*M+1]*Ham_.VL_k[i+M]) * GA[i]) * dFD(w,T);
        
        for (int j = 0; j < M; j++){
            //YX += dw * real(Ham_.VY_RR[i+j*M] * GR[i] * Ham_.VX_LL[j+M*i] * GA[j]) * dFD(w,T);
            //YX_ += dw * real( (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VY_LR[i+j*M] * GR[i] * Ham_.VX_LR[j+M*i] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j]) * dFD(w,T);
            
            XXX0 += dw * Ham_.VX_RR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[j+M*i] * GR[j] * Ham_.VX_LL[j*(M+1)] * GA[j] * dFD(w,T);
            XXX0_ += dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[j+M*i] * GR[j] * Ham_.VX_LR[j*(M+1)] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
            
            XXX1 += dw * Ham_.VX_RR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LL[j+M*i] * GA[j] * dFD(w,T);
            XXX1_ += dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[i+j*M] * GR[i] * GR[i] * Ham_.VXX_LR[j+M*i] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);

            for (int l = 0; l < M; l++){
                if(j==l){

                }
                else{
                    XXX2 += dw * Ham_.VX_RR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LL[l*M+j] * GA[j] * dFD(w,T);
                    XXX2_ += dw * (Ham_.VR_b[j*M]*Ham_.VR_k[j] + Ham_.VR_b[j*M+1]*Ham_.VR_k[j+M]) * Ham_.VX_LR[i+j*M] * GR[i] * GR[i] * Ham_.VX_LR[M*i+l] * GR[l] * Ham_.VX_LR[l*M+j] * (Ham_.VL_b[j*M]*Ham_.VL_k[j] + Ham_.VL_b[j*M+1]*Ham_.VL_k[j+M]) * GA[j] * dFD(w,T);
                }
            }
            
        }
        
    }
    XXX += imag(XXX0+XXX1+XXX2);
    XXX_ += imag(XXX0_+XXX1_+XXX2_);
};

// Calculate linear optical conductivity and photo-galvaic conductivity with Green function
void Opt_transport(parm parm_, double dw, double w,Ham Ham_,Green Green_,double& XX, double& PVYXi, double& div1, double& div2, double& div3){

    Complex VXG[M*M],VYG[M*M],VYXG[M*M],VXGP[M*M],VXGM[M*M],VYGP[M*M],VYGM[M*M],VXGA[M*M],VYGA[M*M],VXGAP[M*M],
        VXGAM[M*M],VYGAP[M*M],VYGAM[M*M],XA[M*M],YA[M*M],XXA[M*M],YXA[M*M],VYXGA[M*M],VXXG[M*M],VXXGA[M*M];
        
    Prod2<M>(Ham_.VX,Green_.GR,VXG);
    Prod2<M>(Ham_.VY,Green_.GR,VYG);

    Prod2<M>(Ham_.VX,Green_.GA,VXGA);
    Prod2<M>(Ham_.VY,Green_.GA,VYGA);

    Prod2<M>(Ham_.VYX,Green_.GR,VYXG);
    Prod2<M>(Ham_.VYX,Green_.GA,VYXGA);

    Prod2<M>(Ham_.VXX,Green_.GR,VXXG);
    Prod2<M>(Ham_.VXX,Green_.GA,VXXGA);

    Prod2<M>(Ham_.VX,Green_.GRp,VXGP);
    Prod2<M>(Ham_.VX,Green_.GRm,VXGM);
    Prod2<M>(Ham_.VX,Green_.GAp,VXGAP);
    Prod2<M>(Ham_.VX,Green_.GAm,VXGAM);
    
    Prod2<M>(Ham_.VY,Green_.GRp,VYGP);
    Prod2<M>(Ham_.VY,Green_.GRm,VYGM);

    
    Prod2<M>(Ham_.VX,Green_.GRmA,XA);
    Prod2<M>(Ham_.VY,Green_.GRmA,YA);
    Prod2<M>(Ham_.VXX,Green_.GRmA,XXA);
    Prod2<M>(Ham_.VYX,Green_.GRmA,YXA);



    //linear term
    Complex PVX1[M*M],PVX2[M*M],PVX3[M*M],PVX4[M*M];

    Prod2<M>(Ham_.VXX,Green_.GRmA,PVX1);
    Prod2<M>(VXGP,XA,PVX2);
    Prod2<M>(XA,VXGAM,PVX3);
    //Prod2<M>(VXGP,VXGA,PVX4);
    Prod2<M>(VXG,VXGA,PVX4);

    Complex PVX = dw * (Trace_NH<M>(PVX1) + Trace_NH<M>(PVX2) + Trace_NH<M>(PVX3)) * FD(w,parm_.T)/parm_.W;
    //Complex PVX_ = dw * (Trace_NH<M>(PVX4)) * dFD(w,parm_.T)/(parm_.W+2.0*I*parm_.delta);


    //nonlinear term
    Complex PVYX1[M*M],PVYX2[M*M],PVYX3[M*M],PVYX4[M*M],PVYX5[M*M],PVYX6[M*M],PVYX7[M*M],PVYX8[M*M];
    Complex PVYX2Q[M*M],PVYX3Q[M*M],PVYX6Q[M*M],PVYX7Q[M*M],PVYX8Q[M*M];

    Prod2<M>(Ham_.VYXX,Green_.GRmA,PVYX1);
    
    Prod3<M>(Ham_.VYX,Green_.GRp,XA,PVYX2);
    Prod2<M>(YXA,VXGAM,PVYX3);

    Prod3<M>(Ham_.VYX,Green_.GRm,XA,PVYX2Q);
    Prod2<M>(YXA,VXGAP,PVYX3Q);

    Prod2<M>(VYG,XXA,PVYX4);
    Prod3<M>(YA,Ham_.VXX,Green_.GA,PVYX5);

    Prod3<M>(VYG,VXGP,XA,PVYX6);
    Prod3<M>(VYGM,XA,VXGAM,PVYX7);
    Prod3<M>(YA,VXGAP,VXGA,PVYX8);

    Prod3<M>(VYG,VXGM,XA,PVYX6Q);
    Prod3<M>(VYGP,XA,VXGAP,PVYX7Q);
    Prod3<M>(YA,VXGAM,VXGA,PVYX8Q);

    Complex PVYX= dw * (Trace_NH<M>(PVYX1) + Trace_NH<M>(PVYX2) + Trace_NH<M>(PVYX2Q) + Trace_NH<M>(PVYX3) + Trace_NH<M>(PVYX3Q) + Trace_NH<M>(PVYX4)
        + Trace_NH<M>(PVYX5) + Trace_NH<M>(PVYX6) + Trace_NH<M>(PVYX6Q) + Trace_NH<M>(PVYX7) + Trace_NH<M>(PVYX7Q) + Trace_NH<M>(PVYX8)+ Trace_NH<M>(PVYX8Q))* FD(w,parm_.T)/(-parm_.W*parm_.W);
    


    //Check of the divergent term
    Complex PVYX1D[M*M],PVYX3D[M*M],PVYX4D[M*M],PVYX6D[M*M];
    Complex PVX1S[M*M],PVX2S[M*M];
    Complex PVYX2S[M*M],PVYX3S[M*M],PVYX4S[M*M],PVYX5S[M*M],PVYX6S[M*M],PVYX7S[M*M],PVYX8S[M*M],PVYX9S[M*M];

    
    //check 1/w term in linear
    Prod2<M>(Ham_.VXX,Green_.GR,PVX1S);
    Prod2<M>(VXG,VXG,PVX2S);
 
    //Check 1/w^2 term in nonlinear
    Prod2<M>(Ham_.VYXX,Green_.GR,PVYX1D);
    Prod2<M>(VYXG,VXG,PVYX3D);
    Prod2<M>(VYG,VXXG,PVYX4D);
    Prod3<M>(VYG,VXG,VXG,PVYX6D);

    //check 1/w term in nonlinear
    Complex VXdG[M*M],VYdG[M*M];
    Prod2<M>(Ham_.VY,Green_.dGR,VYdG);
    Prod2<M>(Ham_.VX,Green_.dGR,VXdG);

    //Fermi sea
    Prod3<M>(Ham_.VYX,Green_.dGR,VXG,PVYX2S);
    Prod2<M>(VYdG,VXXG,PVYX3S);
    Prod3<M>(VYdG,VXG,VXG,PVYX4S);
    Prod3<M>(VYG,VXdG,VXG,PVYX5S);

    //Fermi surface
    Prod2<M>(VYXG,VXGA,PVYX6S);
    Prod3<M>(VYG,VXG,VXGA,PVYX7S);
    Prod2<M>(VYG,VXXGA,PVYX8S);
    Prod3<M>(VYG,VXGA,VXGA,PVYX9S);

    Complex DIV = dw* (Trace_NH<M>(PVYX1D) + 2.0*Trace_NH<M>(PVYX3D)+ Trace_NH<M>(PVYX4D) + 2.0*Trace_NH<M>(PVYX6D)) * FD(w,parm_.T)/(parm_.W*parm_.W);
    Complex DIV2= dw * (Trace_NH<M>(PVX1S)+Trace_NH<M>(PVX2S))*FD(w,parm_.T)/parm_.W;
    Complex DIV3 = dw * ((Trace_NH<M>(PVYX2S)+Trace_NH<M>(PVYX3S)+Trace_NH<M>(PVYX4S)+Trace_NH<M>(PVYX5S))* FD(w,parm_.T) + (Trace_NH<M>(PVYX6S)+Trace_NH<M>(PVYX7S)+Trace_NH<M>(PVYX8S)+Trace_NH<M>(PVYX9S))* dFD(w,parm_.T))/parm_.W;

    div2 += real(DIV2);
    div1 += imag(DIV);
    div3 += imag(DIV3);

    PVYXi += imag(PVYX);
    XX += real(PVX);
};

//For the case of Photo Galvanic effect
void Opt_RTA_transport_BI(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_){
    Complex Drude = 0;
    Complex DrudeL = 0;
    Complex Inj=0;
    Complex BCD =0;
    Complex BC[2];
    for (int i = 0; i < M; i++){
        DrudeL += Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*(M+1)] * (-dFD(Ham_.EN[i],parm_.T)) /(-I*parm_.W + parm_.delta);
        Drude += 4.0* (Ham_.VY_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]* ddFD(Ham_.EN[i],parm_.T) )/(-parm_.W*parm_.W-parm_.delta*parm_.delta);
        for (int j = 0; j < M; j++){
            if(i==j){
                Drude += 4.0* Ham_.VY_LR[i*(M+1)]*Ham_.VXX_LR[i*M+j]* dFD(Ham_.EN[i],parm_.T)/(-parm_.W*parm_.W-parm_.delta*parm_.delta);
            }
            else{
                DrudeL += Ham_.VX_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T)) /(-I*(parm_.W+Ham_.EN[i]-Ham_.EN[j]) + parm_.delta)/(Ham_.EN[i]-Ham_.EN[j]+I* parm_.delta *parm_.W_MAX );
                Drude += 4.0 * Ham_.VY_LR[i*(M+1)]*Ham_.VX_LR[i*M+j]*Ham_.VX_LR[j*M+i]/(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)* dFD(Ham_.EN[i],parm_.T)/(-parm_.W*parm_.W-parm_.delta*parm_.delta);

                Inj += 4.0 * (Ham_.VY_LR[i*(M+1)]-Ham_.VY_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX))*Ham_.VX_LR[i*M+j]*Ham_.VX_LR[j*M+i]
                    *(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);
                
                BC[i] += (Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
            }
        }
        BCD += -2.0 * I * parm_.delta * BC[i] *(Ham_.VX_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
    }

    DrudeL_ += real(DrudeL);
    Drude_ += real(Drude);
    BCD_ += real(BCD);
    Inj_ += real(Inj);

};

//Calculate Linear conductivity, non-reciprocal conductivity, Injection current, and nonlinear Hall effect with RDM methods under RTA 
void Opt_RTA_transport_BI2(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_){
    Complex Drude = 0;
    Complex DrudeL = 0;
    Complex Inj=0;
    Complex BCD =0;
    Complex BC[2];
    for (int i = 0; i < M; i++){
        DrudeL += Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*(M+1)] * (-dFD(Ham_.EN[i],parm_.T)) /(-I*parm_.W + parm_.delta);
        Drude += 4.0* (Ham_.VY_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]* ddFD(Ham_.EN[i],parm_.T) )/(-parm_.W*parm_.W-parm_.delta*parm_.delta);
        for (int j = 0; j < M; j++){
            if(i==j){
                Drude += 4.0* Ham_.VY_LR[i*(M+1)]*Ham_.VXX_LR[i*M+j]* dFD(Ham_.EN[i],parm_.T)/(-parm_.W*parm_.W-parm_.delta*parm_.delta);
            }
            else{
                DrudeL += Ham_.VX_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T)) /(-I*(parm_.W+Ham_.EN[i]-Ham_.EN[j]) + parm_.delta)/(Ham_.EN[i]-Ham_.EN[j]+I* parm_.delta *parm_.W_MAX );
                Drude += 4.0 * Ham_.VY_LR[i*(M+1)]*Ham_.VX_LR[i*M+j]*Ham_.VX_LR[j*M+i]/(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)* dFD(Ham_.EN[i],parm_.T)/(-parm_.W*parm_.W-parm_.delta*parm_.delta);
                Inj += 4.0 * (Ham_.VY_LR[i*(M+1)]-Ham_.VY_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX))*Ham_.VX_LR[i*M+j]*Ham_.VX_LR[j*M+i]
                    *(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);
                
                BC[i] += (Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta/2.0)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta/2.0));
            }
        }
        BCD += -4.0 * I * parm_.delta * BC[i] *(Ham_.VX_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
    }

    DrudeL_ += real(DrudeL);
    Drude_ += real(Drude);
    BCD_ += real(BCD);
    Inj_ += real(Inj);
};

//Calculate Injection current with RDM methods under RTA 
void Opt_RTA_transport_BI3(parm parm_,Ham Ham_,double Inj_[3]){    
    //Complex DrudeL = 0;
    Complex Inj[3]={0,0,0};
    //Complex BCD[3] ={0,0,0};
    //Complex BC[3][2];
    for (int i = 0; i < M; i++){
        //DrudeL += Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*(M+1)] * (-dFD(Ham_.EN[i],parm_.T)) /(-I*parm_.W + parm_.delta);
        
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                //DrudeL += Ham_.VX_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T)) /(-I*(parm_.W+Ham_.EN[i]-Ham_.EN[j]) + parm_.delta)/(Ham_.EN[i]-Ham_.EN[j]+I* parm_.delta *parm_.W_MAX );
                
                Inj[2] += parm_.delta* 4.0 * (Ham_.VZ_LR[i*(M+1)]-Ham_.VZ_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j])*(Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta)*(Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i]-Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);

                Inj[0] += parm_.delta* 4.0 * (Ham_.VX_LR[i*(M+1)]-Ham_.VX_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j])*(Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta)*(Ham_.VY_LR[i*M+j]*Ham_.VZ_LR[j*M+i]-Ham_.VZ_LR[i*M+j]*Ham_.VY_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);

                Inj[1] += parm_.delta* 4.0 * (Ham_.VY_LR[i*(M+1)]-Ham_.VY_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j])*(Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta)*(Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VZ_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);
                
                //BC[2][i] += (Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
                //BC[0][i] += (Ham_.VZ_LR[i*M+j]*Ham_.VY_LR[j*M+i]-Ham_.VY_LR[i*M+j]*Ham_.VZ_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
                //BC[1][i] += (Ham_.VX_LR[i*M+j]*Ham_.VZ_LR[j*M+i]-Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
            }
        }
        //BCD[0] += -2.0 * I * parm_.delta * BC[0][i] *(Ham_.VX_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
        //BCD[1] += -2.0 * I * parm_.delta * BC[1][i] *(Ham_.VY_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
        //BCD[2] += -2.0 * I * parm_.delta * BC[2][i] *(Ham_.VZ_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
    }

    //DrudeL_ += real(DrudeL);

    for (int i = 0; i < 3; i++){
        Inj_[i] += imag(Inj[i]);
        //BCD_[i] += real(BCD[i]);
    }

};

//Calculate Linear conductivity, Injection current, and nonlinear Hall effect with RDM methods under RTA 
void Opt_RTA_transport_BI4(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3]){
    
    Complex DrudeL = 0;
    Complex Inj[3]={0,0,0};
    Complex BCD[3] ={0,0,0};
    Complex BC[3][2];
    for (int i = 0; i < M; i++){
        DrudeL += Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*(M+1)] * (-dFD(Ham_.EN[i],parm_.T)) /(-I*parm_.W + parm_.delta);
        
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                DrudeL += Ham_.VX_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T)) /(-I*(parm_.W+Ham_.EN[i]-Ham_.EN[j]) + parm_.delta)/(Ham_.EN[i]-Ham_.EN[j]+I* parm_.delta *parm_.W_MAX );
                
                if(fabs(parm_.W + Ham_.EN[i]-Ham_.EN[j])<parm_.delta){
                    Inj[2] += parm_.delta* 4.0 * (Ham_.VZ_LR[i*(M+1)]-Ham_.VZ_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta))*(Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i]-Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T));

                    Inj[0] += parm_.delta* 4.0 * (Ham_.VX_LR[i*(M+1)]-Ham_.VX_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta))*(Ham_.VY_LR[i*M+j]*Ham_.VZ_LR[j*M+i]-Ham_.VZ_LR[i*M+j]*Ham_.VY_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T));

                    Inj[1] += parm_.delta* 4.0 * (Ham_.VY_LR[i*(M+1)]-Ham_.VY_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta))*(Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VZ_LR[j*M+i])*(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T));
                }
                
                BC[2][i] += (Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
                BC[0][i] += (Ham_.VZ_LR[i*M+j]*Ham_.VY_LR[j*M+i]-Ham_.VY_LR[i*M+j]*Ham_.VZ_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
                BC[1][i] += (Ham_.VX_LR[i*M+j]*Ham_.VZ_LR[j*M+i]-Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX));
            }
        }
        BCD[0] += -2.0 * I * parm_.delta * BC[0][i] *(Ham_.VX_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
        BCD[1] += -2.0 * I * parm_.delta * BC[1][i] *(Ham_.VY_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
        BCD[2] += -2.0 * I * parm_.delta * BC[2][i] *(Ham_.VZ_LR[i*(M+1)] * dFD(Ham_.EN[i],parm_.T))/ (-parm_.W*parm_.W-parm_.delta*parm_.delta);
    }

    DrudeL_ += real(DrudeL);

    for (int i = 0; i < 3; i++){
        Inj_[i] += imag(Inj[i]);
        BCD_[i] += real(BCD[i]);
    }

};


//Calculate Linear conductivity, Injection current, and nonlinear Hall effect by Green function methods without non-Hermitian effect
void Opt_Green_transport_BI(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3] ,double Inj2_[3]){
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
        //Complex Drude = 0;
        Complex DrudeL = 0;
        Complex Inj[3]={0,0,0};
        Complex Inj2[3]={0,0,0};
        Complex BCD[3] ={0,0,0};
        for (int i = 0; i < M; i++){
            DrudeL += Ham_.VX_LR[i*(M+1)]/((w -Ham_.EN[i])*(w -Ham_.EN[i])- (parm_.W+I*parm_.delta)*(parm_.W+I*parm_.delta))* Ham_.VX_LR[i*(M+1)]*(1.0/(w -Ham_.EN[i]+I*parm_.delta)-1.0/(w -Ham_.EN[i]-I*parm_.delta))* FD(w,parm_.T);

            for (int j = 0; j < M; j++){
                if(i==j){
                }
                else{                    

                    Complex A0 = I*parm_.delta*(FD(w+parm_.W,parm_.T)-FD(w,parm_.T))/((w-Ham_.EN[i])*(w-Ham_.EN[i])+parm_.delta*parm_.delta)/((w+parm_.W-Ham_.EN[j])*(w+parm_.W-Ham_.EN[j])+parm_.delta*parm_.delta)/(-parm_.W*parm_.W);

                    //Complex A1 = FD(w,parm_.T)/((w-Ham_.EN[i]+I*parm_.delta)*(w-Ham_.EN[i]+I*parm_.delta))/(w+parm_.W-Ham_.EN[j]+I*parm_.delta)/(-parm_.W*parm_.W);

                    Complex A1 = 0; Complex B1 = 0;
                    Complex B0 = I*parm_.delta*(FD(w-parm_.W,parm_.T)-FD(w,parm_.T))/((w-Ham_.EN[i])*(w-Ham_.EN[i])+parm_.delta*parm_.delta)/((w-parm_.W-Ham_.EN[j])*(w-parm_.W-Ham_.EN[j])+parm_.delta*parm_.delta)/(-parm_.W*parm_.W);

                    //Complex B1 = FD(w,parm_.T)/((w-Ham_.EN[i]+I*parm_.delta)*(w-Ham_.EN[i]+I*parm_.delta))/(w-parm_.W-Ham_.EN[j]+I*parm_.delta)/(-parm_.W*parm_.W);

                    Inj[2] += 2.0 * Ham_.VZ_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] * (A0 + A1) + 2.0 * Ham_.VZ_LR[i*(M+1)] * Ham_.VY_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (B0 + B1);

                    Inj[0] += 2.0 * Ham_.VX_LR[i*(M+1)] * Ham_.VY_LR[i*M+j] * Ham_.VZ_LR[j*M+i] * (A0 + A1) + 2.0 * Ham_.VX_LR[i*(M+1)] * Ham_.VZ_LR[i*M+j] * Ham_.VY_LR[j*M+i] * (B0 + B1);

                    Inj[1] += 2.0 * Ham_.VY_LR[i*(M+1)] * Ham_.VZ_LR[i*M+j] * Ham_.VX_LR[j*M+i] * (A0 + A1) + 2.0 * Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VZ_LR[j*M+i] * (B0 + B1);
                }
            }
        }

        //DrudeL_ += dw * real(DrudeL)/(2.0*pi);
        for (int i = 0; i < 3; i++){
            Inj_[i] += dw * real(Inj[i])/(2.0*pi);
            //Inj2_[i] += dw * imag(Inj2[i])/(2.0*pi);
            //BCD_[i] += dw * imag(BCD[i])/(2.0*pi);
        }
    } 
};

void Opt_Green_transport_BI_div(parm parm_,Ham Ham_,double w,double dk2, double div_[3]){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            div_[0] += 2.0 * dk2* imag(Ham_.VX_LR[i*(M+1)]/(w-Ham_.EN[i]+I*parm_.delta)*(Ham_.VY_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VZ_LR[j*M+i]+Ham_.VZ_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VY_LR[j*M+i])/(w-Ham_.EN[i]+I*parm_.delta))/(4.0*pi*pi*pi);

            div_[1] += 2.0 * dk2* imag(Ham_.VY_LR[i*(M+1)]/(w-Ham_.EN[i]+I*parm_.delta)*(Ham_.VZ_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VX_LR[j*M+i]+Ham_.VX_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VZ_LR[j*M+i])/(w-Ham_.EN[i]+I*parm_.delta))/(4.0*pi*pi*pi);

            div_[2] += 2.0 * dk2* imag(Ham_.VZ_LR[i*(M+1)]/(w-Ham_.EN[i]+I*parm_.delta)*(Ham_.VX_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VY_LR[j*M+i]+Ham_.VY_LR[i*M+j]/(w-Ham_.EN[j]+I*parm_.delta)*Ham_.VX_LR[j*M+i])/(w-Ham_.EN[i]+I*parm_.delta))/(4.0*pi*pi*pi);
        }
    }
};


//Calculate intrinsic Fermi surface term with RDM methods under RTA
void Opt_IFS_BI(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                double Del = Ham_.EN[i]-Ham_.EN[j];
                Complex Gm = (Ham_.VX_LR[i+M*j]*Ham_.VY_LR[j+M*i] + Ham_.VY_LR[i+M*j]*Ham_.VX_LR[j+M*i])/(Del*Del + parm_.delta*parm_.delta);
                IFS_xxy += real(Gm/((parm_.W - Del)*(parm_.W - Del) + parm_.delta*parm_.delta)*Del*Ham_.VX_LR[i*(M+1)]*dFD(Ham_.EN[i],parm_.T));
                IFS_yxy += real(Gm/((parm_.W - Del)*(parm_.W - Del) + parm_.delta*parm_.delta)*Del*Ham_.VY_LR[i*(M+1)]*dFD(Ham_.EN[i],parm_.T));
            }
        }
        
    }
    
}

//Calculate intrinsic Fermi surface term by the Green function methods without non-Hermitian effect
void Opt_IFS_Green(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy){
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        Complex IFS1_xxy=0;
        Complex IFS2_xxy=0;
        Complex IFS1_yxy=0;
        Complex IFS2_yxy=0;
        for (int i = 0; i < M; i++){
            for (int j = 0; j < M; j++){
                if(i==j){
                }
                else{
                    double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
                    Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
                    Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);
                    Complex GjR = 1.0/(w-Ham_.EN[j]+I*parm_.delta/2.0);
                    Complex GjRm = 1.0/(w - parm_.W - Ham_.EN[j] + I*parm_.delta/2.0);
                    IFS1_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VX_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                    IFS1_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VY_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                    IFS2_xxy += FD(w,parm_.T) * (Ham_.VXX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYX_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );

                    IFS2_yxy += FD(w,parm_.T) * (Ham_.VYX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYY_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );
                    
                }
            }
            
        }
        IFS_xxy += -dw * imag(IFS1_xxy) / (2*pi*(parm_.W*parm_.W+0.25*parm_.delta*parm_.delta));
        IFS_yxy += -dw * imag(IFS1_yxy) / (2*pi*(parm_.W*parm_.W+0.25*parm_.delta*parm_.delta));
    }
    
}

void Opt_IFS_Green_w(parm parm_, Ham Ham_,double w, double& IFS_xxy, double& IFS_yxy){
    
    //double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    //for (int ww = 0; ww < parm_.W_SIZE; ww++){
    Complex IFS1_xxy=0;
    Complex IFS2_xxy=0;
    Complex IFS1_yxy=0;
    Complex IFS2_yxy=0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                //double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
                Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
                Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);
                Complex GjR = 1.0/(w-Ham_.EN[j]+I*parm_.delta/2.0);
                Complex GjRm = 1.0/(w - parm_.W - Ham_.EN[j] + I*parm_.delta/2.0);
                IFS1_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                    *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                            + Ham_.VX_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                    * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                IFS1_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                    *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                            + Ham_.VY_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                    * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                IFS2_xxy += FD(w,parm_.T) * (Ham_.VXX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                            + Ham_.VYX_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );

                IFS2_yxy += FD(w,parm_.T) * (Ham_.VYX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                            + Ham_.VYY_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );
                
            }
        }
        
    }
    IFS_xxy += -imag(IFS1_xxy+IFS2_xxy) / (2*pi*(parm_.W*parm_.W));
    IFS_yxy += -imag(IFS1_yxy+IFS2_yxy) / (2*pi*(parm_.W*parm_.W));
    //}
    
}


void Opt_IFS_Green_w_div(parm parm_, Ham Ham_,double w, double& div_xxy, double& div_yxy){
    
    //double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    //for (int ww = 0; ww < parm_.W_SIZE; ww++){
    Complex div0_xxy=0;
    Complex div0_yxy=0;
    for (int i = 0; i < M; i++){

        Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
        Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);

        div0_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VYX_LR[i*M+i] *( (GiR-conj(GiR))*(GiR+conj(GiR)) )
                                    + Ham_.VYXX_LR[i*M+i]  * (GiR-conj(GiR)) );

        div0_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VYX_LR[i*M+i] *( (GiR-conj(GiR))*(GiR+conj(GiR)) )
                                    + Ham_.VYYX_LR[i*M+i]  * (GiR-conj(GiR)) );
        
    }
    div_xxy += -imag(div0_xxy) / (2*pi*(parm_.W*parm_.W));
    div_yxy += -imag(div0_yxy) / (2*pi*(parm_.W*parm_.W));
    //}
    
}

void Opt_IFS_Green_test(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy){
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        Complex IFS1_xxy=0;
        Complex IFS2_xxy=0;
        Complex IFS1_yxy=0;
        Complex IFS2_yxy=0;
        for (int i = 0; i < M; i++){
            for (int j = 0; j < M; j++){
                if(i==j){
                }
                else{
                    double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
                    double del = Ham_.EN[i] - Ham_.EN[j];
                    Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
                    Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);
                    Complex GjR = 1.0/(w-Ham_.EN[j]+I*parm_.delta/2.0);
                    Complex GjRm = 1.0/(w - parm_.W - Ham_.EN[j] + I*parm_.delta/2.0);
                    IFS1_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VX_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) )/(del*del+parm_.delta*parm_.delta/4.0);

                    IFS1_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VY_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) )/(del*del+parm_.delta*parm_.delta/4.0);

                    IFS2_xxy += FD(w,parm_.T) * (Ham_.VXX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYX_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) )/(del*del+parm_.delta*parm_.delta/4.0);

                    IFS2_yxy += FD(w,parm_.T) * (Ham_.VYX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYY_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) )/(del*del+parm_.delta*parm_.delta/4.0);
                    
                }
            }
            
        }
        IFS_xxy += -dw * imag(IFS1_xxy+IFS2_xxy) / (2*pi);
        IFS_yxy += -dw * imag(IFS1_yxy+IFS2_yxy) / (2*pi);
    }
    
}

void Opt_IFS_Green_test2(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy){
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        Complex IFS1_xxy=0;
        Complex IFS2_xxy=0;
        Complex IFS1_yxy=0;
        Complex IFS2_yxy=0;
        for (int i = 0; i < M; i++){
            for (int j = 0; j < M; j++){
                if(i==j){
                }
                else{
                    double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
                    double del = Ham_.EN[i] - Ham_.EN[j];
                    Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
                    Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);
                    Complex GjR = 1.0/(w-Ham_.EN[j]+I*parm_.delta/2.0);
                    Complex GjRm = 1.0/(w - parm_.W - Ham_.EN[j] + I*parm_.delta/2.0);
                    IFS1_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VX_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) )/(del*del+parm_.delta*parm_.delta);

                    IFS1_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                     *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                                + Ham_.VY_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                     * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) )/(del*del+parm_.delta*parm_.delta);

                    IFS2_xxy += FD(w,parm_.T) * (Ham_.VXX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYX_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) )/(del*del+parm_.delta*parm_.delta);

                    IFS2_yxy += FD(w,parm_.T) * (Ham_.VYX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                                + Ham_.VYY_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) )/(del*del+parm_.delta*parm_.delta);
                    
                }
            }
            
        }
        IFS_xxy += -dw * imag(IFS1_xxy+IFS2_xxy) / (2*pi);
        IFS_yxy += -dw * imag(IFS1_yxy+IFS2_yxy) / (2*pi);
    }
    
}

/*
void Opt_IFS_NH_BI(parm parm_, Ham Ham_, double& IFS_RA, double& IFS_RR, double& div1, double& div2){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                double Del = Ham_.EN[i]-Ham_EN[j];
            }
        }
        
    }
    
}*/

void Opt_Gyration_BI(parm parm_, Ham Ham_, double& gyro_xxy, double& gyro_yxy){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                double e_ij = Ham_.EN[i]-Ham_.EN[j];
                Complex delvx_ij = Ham_.VX_LR[i*(M+1)] - Ham_.VX_LR[j*(M+1)];
                Complex delvy_ij = Ham_.VY_LR[i*(M+1)] - Ham_.VY_LR[j*(M+1)];
                double f_ij = FD(Ham_.EN[i],parm_.T) - FD(Ham_.EN[j],parm_.T);
                double del = parm_.delta/((parm_.W-e_ij)*(parm_.W-e_ij) + parm_.delta*parm_.delta);
                for (int l = 0; l < M; l++){
                    double e_il = Ham_.EN[i]-Ham_.EN[l];
                    double e_lj = Ham_.EN[l]-Ham_.EN[j];
                    gyro_xxy += pi * real((-delvx_ij*Ham_.VX_LR[M*i+j]/(e_ij*e_ij+parm_.delta*parm_.delta) 
                                             + e_ij*(Ham_.VXX_LR[M*i+j]
                                                 + Ham_.VX_LR[M*i+l]*e_il*Ham_.VX_LR[M*l+j]/(e_il*e_il + parm_.delta*parm_.delta)
                                                 -Ham_.VX_LR[M*i+l]*e_lj*Ham_.VX_LR[M*l+j]/(e_lj*e_lj + parm_.delta*parm_.delta) )/(e_ij*e_ij + parm_.delta*parm_.delta))
                                            *Ham_.VY_LR[j*M+i]*(-e_ij)/(e_ij*e_ij + parm_.delta*parm_.delta)  ) *f_ij*del;

                    gyro_xxy += pi * real( (delvx_ij*Ham_.VY_LR[M*i+j]/(e_ij*e_ij+parm_.delta*parm_.delta)
                                            - e_ij*(Ham_.VYX_LR[M*i+j]
                                                 + Ham_.VX_LR[M*i+l]*e_il*Ham_.VY_LR[M*l+j]/(e_il*e_il + parm_.delta*parm_.delta)
                                                 -Ham_.VX_LR[M*i+l]*e_lj*Ham_.VY_LR[M*l+j]/(e_lj*e_lj + parm_.delta*parm_.delta) )/(e_ij*e_ij + parm_.delta*parm_.delta))
                                            *Ham_.VX_LR[j*M+i]*(-e_ij)/(e_ij*e_ij + parm_.delta*parm_.delta) )*f_ij*del;

                    gyro_yxy += pi * real((-delvy_ij*Ham_.VX_LR[M*i+j]/(e_ij*e_ij+parm_.delta*parm_.delta) 
                                             + e_ij*(Ham_.VYX_LR[M*i+j]
                                                 + Ham_.VY_LR[M*i+l]*e_il*Ham_.VX_LR[M*l+j]/(e_il*e_il + parm_.delta*parm_.delta)
                                                 -Ham_.VY_LR[M*i+l]*e_lj*Ham_.VX_LR[M*l+j]/(e_lj*e_lj + parm_.delta*parm_.delta) )/(e_ij*e_ij + parm_.delta*parm_.delta))
                                            *Ham_.VY_LR[j*M+i]*(-e_ij)/(e_ij*e_ij + parm_.delta*parm_.delta)  ) *f_ij*del;

                    gyro_yxy += pi * real( (delvy_ij*Ham_.VY_LR[M*i+j]/(e_ij*e_ij+parm_.delta*parm_.delta)
                                            - e_ij*(Ham_.VYY_LR[M*i+j]
                                                 + Ham_.VY_LR[M*i+l]*e_il*Ham_.VY_LR[M*l+j]/(e_il*e_il + parm_.delta*parm_.delta)
                                                 -Ham_.VY_LR[M*i+l]*e_lj*Ham_.VY_LR[M*l+j]/(e_lj*e_lj + parm_.delta*parm_.delta) )/(e_ij*e_ij + parm_.delta*parm_.delta))
                                            *Ham_.VX_LR[j*M+i]*(-e_ij)/(e_ij*e_ij + parm_.delta*parm_.delta) )*f_ij*del;
                }
                
            }
        }
    }
}


void Opt_Gyration_BI2(parm parm_, Ham Ham_, double& gyro_xxy){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                double e_ij = Ham_.EN[i]-Ham_.EN[j];
                Complex delvx_ij = Ham_.VX_LR[i*(M+1)] - Ham_.VX_LR[j*(M+1)];
                Complex delvy_ij = Ham_.VY_LR[i*(M+1)] - Ham_.VY_LR[j*(M+1)];
                double f_ij = FD(Ham_.EN[i],parm_.T) - FD(Ham_.EN[j],parm_.T);
                double del = parm_.delta/((parm_.W-e_ij)*(parm_.W-e_ij) + parm_.delta*parm_.delta);
                double delm = parm_.delta/((-parm_.W-e_ij)*(-parm_.W-e_ij) + parm_.delta*parm_.delta);
                gyro_xxy += real(Ham_.VXX_LR[i*M+j]*Ham_.VY_LR[j*M+i]/(e_ij*e_ij + parm_.delta*parm_.delta)) *f_ij*del;

                gyro_xxy += real(Ham_.VYX_LR[i*M+j]*Ham_.VX_LR[j*M+i]/(e_ij*e_ij + parm_.delta*parm_.delta)) *f_ij*delm;
                
            }
        }
    }
}

void Opt_Gyration_BI3(parm parm_, Ham Ham_, double& gyro_xxy, double& gyro_xxy_im){
    gyro_xxy = 0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            //if(i==j){
            //}
            //else{
                double e_ij = Ham_.EN[i]-Ham_.EN[j];
                //Complex delvx_ij = Ham_.VX_LR[i*(M+1)] - Ham_.VX_LR[j*(M+1)];
                //Complex delvy_ij = Ham_.VY_LR[i*(M+1)] - Ham_.VY_LR[j*(M+1)];
                double f_ij = FD(Ham_.EN[i],parm_.T) - FD(Ham_.EN[j],parm_.T);
                double del = parm_.delta/((parm_.W+e_ij)*(parm_.W+e_ij) + parm_.delta*parm_.delta);
                double delm = parm_.delta/((-parm_.W+e_ij)*(-parm_.W+e_ij) + parm_.delta*parm_.delta);
                Complex deli = 1.0/((parm_.W+e_ij) + I*parm_.delta);
                Complex delmi = 1.0/((-parm_.W+e_ij) + I*parm_.delta);
                gyro_xxy += real(Ham_.VXX_LR[i*M+j]*Ham_.VY_LR[j*M+i]/(parm_.W*parm_.W)) *f_ij*del;

                gyro_xxy += real(Ham_.VYX_LR[i*M+j]*Ham_.VX_LR[j*M+i]/(parm_.W*parm_.W)) *f_ij*delm;
                gyro_xxy_im += real(Ham_.VXX_LR[i*M+j]*Ham_.VY_LR[j*M+i]/(parm_.W*parm_.W)*deli) *f_ij;

                gyro_xxy_im += real(Ham_.VYX_LR[i*M+j]*Ham_.VX_LR[j*M+i]/(parm_.W*parm_.W)*delmi) *f_ij;
                
            //}
        }
    }
}

void Opt_Gyration_Green(parm parm_, Ham Ham_, double w, double& gyro_xxy, double& gyro_yxy){
    //double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    //for (int ww = 0; ww < parm_.W_SIZE; ww++){
    Complex gyro1_xxy=0;
    Complex gyro2_xxy=0;
    Complex gyro1_yxy=0;
    Complex gyro2_yxy=0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            if(i==j){
            }
            else{
                //double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
                Complex GiR = 1.0/(w-Ham_.EN[i]+I*parm_.delta/2.0);
                Complex GiRp = 1.0/(w + parm_.W - Ham_.EN[i] + I*parm_.delta/2.0);
                Complex GjR = 1.0/(w-Ham_.EN[j]+I*parm_.delta/2.0);
                Complex GjRm = 1.0/(w - parm_.W - Ham_.EN[j] + I*parm_.delta/2.0);
                gyro1_xxy += FD(w,parm_.T) * (Ham_.VX_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                    *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                            - Ham_.VX_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                    * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                gyro1_yxy += FD(w,parm_.T) * (Ham_.VY_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i]
                                                    *( GiR*GjRm*(GiR-conj(GiR)) + GiRp*(GjR-conj(GjR))*conj(GiRp) + conj(GiR)*conj(GjRm)*(GiR-conj(GiR)))
                                            - Ham_.VY_LR[j*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] 
                                                    * ( GjR*GiRp*(GjR-conj(GjR)) + GjRm*(GiR-conj(GiR))*conj(GjRm) + conj(GjR)*conj(GiRp)*(GjR-conj(GjR)) ) );

                gyro2_xxy += FD(w,parm_.T) * (Ham_.VXX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                            - Ham_.VYX_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );

                gyro2_yxy += FD(w,parm_.T) * (Ham_.VYX_LR[i*M+j] * Ham_.VY_LR[j*M+i] *( GjRm*(GiR-conj(GiR)) + (GjR-conj(GjR))*conj(GiRp) )
                                            - Ham_.VYY_LR[j*M+i] * Ham_.VX_LR[i*M+j]  * ( GiRp*(GjR-conj(GjR)) + (GiR-conj(GiR))*conj(GjRm) ) );
                
            }
        }
        
    }
    gyro_xxy += -real(gyro1_xxy+gyro2_xxy) / (2*pi*(parm_.W*parm_.W));
    gyro_yxy += -real(gyro1_yxy+gyro2_yxy) / (2*pi*(parm_.W*parm_.W));
    //}
    
}


// Calculate photo-galvaic conductivity under circular polarized light with Green function
void Opt_Circular_Green(parm parm_,Ham Ham_,Green Green_,double w, double& PVCP1, double& PVCP2, double& PVCPWW, double& PVCP){

    double ff = FD(w,parm_.T);
    Complex VXXGRP[M*M], VYGRA[M*M], VXXGRA[M*M], VYGAM[M*M], VXYGRM[M*M], VXGRA[M*M], VXYGRA[M*M], VXGAP[M*M];
    Prod2<M>(Ham_.VXX,Green_.GRp,VXXGRP);
    Prod2<M>(Ham_.VY, Green_.GRmA,VYGRA);
    Prod2<M>(Ham_.VXX,Green_.GRmA,VXXGRA);
    Prod2<M>(Ham_.VY, Green_.GAm,VYGAM);

    Prod2<M>(Ham_.VYX,Green_.GRm,VXYGRM);
    Prod2<M>(Ham_.VX, Green_.GRmA,VXGRA);
    Prod2<M>(Ham_.VYX,Green_.GRmA,VXYGRA);
    Prod2<M>(Ham_.VX, Green_.GAp,VXGAP);

    Complex PV3_1[M*M], PV4_1[M*M],PV3_2[M*M], PV4_2[M*M];
    Prod2<M>(VXXGRP, VYGRA, PV3_1);
    Prod2<M>(VXXGRA, VYGAM, PV3_2);
    Prod2<M>(VXYGRM, VXGRA, PV4_1);
    Prod2<M>(VXYGRA, VXGAP, PV4_2);

    Complex VXGR[M*M], VXGRM[M*M], VYGA[M*M], VXGRP[M*M];
    Prod2<M>(Ham_.VX, Green_.GR, VXGR);
    Prod2<M>(Ham_.VX, Green_.GRp, VXGRP);
    Prod2<M>(Ham_.VX, Green_.GRm, VXGRM);
    Prod2<M>(Ham_.VY, Green_.GA, VYGA);

    Complex PV5_1[M*M], PV5_2[M*M], PV5_3[M*M];
    Prod3<M>(VXGR, VXGRP, VYGRA, PV5_1);
    Prod3<M>(VXGRM,VXGRA,VYGAM, PV5_2);
    Prod3<M>(VXGRA, VXGAP, VYGA, PV5_3);
    double ddk = 4.0 * pi * pi /(parm_.K_SIZE*parm_.K_SIZE);
    //double ddw = 2.0 * parm_.W_MAX/parm_.W_SIZE; 
    double pv1 = ff*ddk*(Trace_H<M>(PV3_1)+Trace_H<M>(PV3_2)+Trace_H<M>(PV4_1)+Trace_H<M>(PV4_2));
    double pv2 = ff*ddk*(Trace_H<M>(PV5_1)+Trace_H<M>(PV5_2)+Trace_H<M>(PV5_3));
    PVCP1 += pv1;
    PVCP2 += pv2;
    PVCPWW += pv1+pv2;
    PVCP += (pv1+pv2)/(parm_.W*parm_.W);

};

void Opt_Circular_Green_BI(parm parm_,Ham Ham_,double w, double& PVCP1, double& PVCP){

    double PV1 = 0;
    double ddk = 4.0 * pi * pi /(parm_.K_SIZE*parm_.K_SIZE);
    double ff = FD(w,parm_.T);
    for (int i = 0; i < M; i++){
        Complex GRPi = 1.0/(w+parm_.W-Ham_.EN[i]+I*parm_.delta);
        Complex GRMi = 1.0/(w-parm_.W-Ham_.EN[i]+I*parm_.delta);
        Complex GRmAi = 1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta);
        for (int j = 0; j < M; j++){
            Complex GAPj = 1.0/(w+parm_.W-Ham_.EN[j]-I*parm_.delta);
            Complex GAMj = 1.0/(w-parm_.W-Ham_.EN[j]-I*parm_.delta);
            Complex GRmAj = 1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta);
            PV1 += real(Ham_.VXX_LR[j*M+i]*GRPi*Ham_.VY_LR[i*M+j]*GRmAj);
            PV1 += real(Ham_.VXX_LR[j*M+i]*GRmAi*Ham_.VY_LR[i*M+j]*GAMj);
            PV1 += real(Ham_.VYX_LR[j*M+i]*GRMi*Ham_.VX_LR[i*M+j]*GRmAj);
            PV1 += real(Ham_.VYX_LR[j*M+i]*GRmAi*Ham_.VX_LR[i*M+j]*GAPj);
        }
    }
    PVCP1 += ddk*ff*PV1;
    PVCP += ddk*ff*PV1/(parm_.W*parm_.W);
};
