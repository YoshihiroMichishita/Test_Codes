#include "transport.hpp"

using namespace std;

double FD(double w, double T){
    double f = 1.0/(1.0+exp(w/T));
    return f;
};

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

double Spectral(parm parm_,double w, double k[D], double re[Mf], double im[Mf]){
    Complex G[M*M],H[M*M];
    H_mom(parm_,k,H);
    GreenR_mom(parm_,w,im,re,H,G);
    double A = -imag(Trace_NH<M>(G))/pi;
    return A;
};

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

void Linear_transport_NH(double dw, double w,double T,Ham Ham_, double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi,
 double& YXXi2, double& XX_, double& YX_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_){
    
    Complex GR[M],GA[M];
    for (int i = 0; i < M; i++){
        GR[i] = 1.0/(w - Ham_.E_NH[i]);
        GA[i] = conj(GR[i]);
    }
    Complex YXX0=0,YXX0_=0,YXX1=0,YXX1_=0;

    for (int i = 0; i < M; i++){
        XX += dw * real(Ham_.VX_RR[i*(M+1)] * GR[i] * Ham_.VX_LL[i*(M+1)] * GA[i]) * dFD(w,T);
        XX_ += dw * real( (Ham_.VR_b[M*i]*Ham_.VR_k[i] + Ham_.VR_b[M*i+1]*Ham_.VR_k[i+M]) * Ham_.VX_LR[i*(M+1)] * GR[i] * Ham_.VX_LR[i*(M+1)] * (Ham_.VL_b[M*i]*Ham_.VL_k[i] + Ham_.VL_b[M*i+1]*Ham_.VL_k[i+M]) * GA[i]) * dFD(w,T);

        for (int j = 0; j < M; j++){
            for (int l = 0; l < M; l++){
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
    YXXi += imag(YXX0);
    YXXr += real(YXX0);
    YXXi_ += imag(YXX0_);
    YXXr_ += real(YXX0_);
    YXXi2 += imag(YXX1);
    YXXr2 += real(YXX1);
    YXXi2_ += imag(YXX1_);
    YXXr2_ += real(YXX1_);
};

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

void Opt_RTA_transport_BI3(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3]){
    
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
                
                Inj[2] += 4.0 * (Ham_.VZ_LR[i*(M+1)]-Ham_.VZ_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX))*Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i]
                    *(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);

                Inj[0] += 4.0 * (Ham_.VX_LR[i*(M+1)]-Ham_.VX_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX))*Ham_.VY_LR[i*M+j]*Ham_.VZ_LR[j*M+i]
                    *(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);

                Inj[1] += 4.0 * (Ham_.VY_LR[i*(M+1)]-Ham_.VY_LR[j*(M+1)])/((Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta)*(Ham_.EN[i]-Ham_.EN[j]+I*parm_.delta*parm_.W_MAX))*Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i]
                    *(FD(Ham_.EN[i],parm_.T)-FD(Ham_.EN[j],parm_.T))/((parm_.W+Ham_.EN[i]-Ham_.EN[j])*(parm_.W+Ham_.EN[i]-Ham_.EN[j])+parm_.delta*parm_.delta);
                
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
        Inj_[i] += real(Inj[i]);
        BCD_[i] += real(BCD[i]);
    }

};

void Opt_Green_transport_BI(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3] ,double Inj2_[3]){
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        double w = parm_.W_MAX * (ww -parm_.W_SIZE/2) * 2.0/parm_.W_SIZE;
        Complex Drude = 0;
        Complex DrudeL = 0;
        Complex Inj[3]={0,0,0};
        Complex Inj2[3]={0,0,0};
        Complex BCD[3] ={0,0,0};
        for (int i = 0; i < M; i++){
            DrudeL += Ham_.VX_LR[i*(M+1)]/((w -Ham_.EN[i])*(w -Ham_.EN[i])- (parm_.W+I*parm_.delta)*(parm_.W+I*parm_.delta))* Ham_.VX_LR[i*(M+1)]*(1.0/(w -Ham_.EN[i]+I*parm_.delta)-1.0/(w -Ham_.EN[i]-I*parm_.delta))* FD(w,parm_.T);

            //Drude += 4.0* (Ham_.VY_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]*Ham_.VX_LR[i*(M+1)]* ddFD(Ham_.EN[i],parm_.T) )/(-parm_.W*parm_.W-parm_.delta*parm_.delta);

            for (int j = 0; j < M; j++){
                if(i==j){
                }
                else{                    
                    Inj[2] += 4.0 * Ham_.VZ_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] * FD(w,parm_.T) *(2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) /(w-Ham_.EN[i]+I*parm_.delta) / (w + parm_.W -Ham_.EN[j] + I*parm_.delta) + (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta)) /(w-parm_.W-Ham_.EN[i]+I*parm_.delta) / (w - parm_.W -Ham_.EN[j] - I*parm_.delta) )/(-parm_.W*parm_.W);

                    Inj[1] += 4.0 * Ham_.VY_LR[i*(M+1)] * Ham_.VZ_LR[i*M+j] * Ham_.VX_LR[j*M+i] * FD(w,parm_.T) *(2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) /(w-Ham_.EN[i]+I*parm_.delta) / (w + parm_.W -Ham_.EN[j] + I*parm_.delta) + (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta)) /(w-parm_.W-Ham_.EN[i]+I*parm_.delta) / (w - parm_.W -Ham_.EN[j] - I*parm_.delta) )/(-parm_.W*parm_.W);

                    Inj[0] += 4.0 * Ham_.VX_LR[i*(M+1)] * Ham_.VY_LR[i*M+j] * Ham_.VZ_LR[j*M+i] * FD(w,parm_.T) *(2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) /(w-Ham_.EN[i]+I*parm_.delta) / (w + parm_.W -Ham_.EN[j] + I*parm_.delta) + (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta)) /(w-parm_.W-Ham_.EN[i]+I*parm_.delta) / (w - parm_.W -Ham_.EN[j] - I*parm_.delta) )/(-parm_.W*parm_.W);

                    Inj2[2] += 2.0 * Ham_.VZ_LR[i*(M+1)] * Ham_.VX_LR[i*M+j] * Ham_.VY_LR[j*M+i] * FD(w,parm_.T) *2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) / (w-Ham_.EN[i]+I*parm_.delta)/(Ham_.EN[j]-Ham_.EN[i]+I*parm_.delta) / (-parm_.W*parm_.W);

                    Inj2[0] += 2.0 * Ham_.VY_LR[i*(M+1)] * Ham_.VZ_LR[i*M+j] * Ham_.VX_LR[j*M+i] * FD(w,parm_.T) *2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) / (w-Ham_.EN[i]+I*parm_.delta)/(Ham_.EN[j]-Ham_.EN[i]+I*parm_.delta) / (-parm_.W*parm_.W);

                    Inj2[1] += 2.0 * Ham_.VX_LR[i*(M+1)] * Ham_.VY_LR[i*M+j] * Ham_.VZ_LR[j*M+i] * FD(w,parm_.T) *2.0*(1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta)) / (w-Ham_.EN[i]+I*parm_.delta)/(Ham_.EN[j]-Ham_.EN[i]+I*parm_.delta) / (-parm_.W*parm_.W);

                    BCD[0] += (Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i])*Ham_.VX_LR[i*(M+1)] * FD(w,parm_.T) *((1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta))*(1.0/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[j] + I*parm_.delta) + 1.0/(w-parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-parm_.W-Ham_.EN[j] + I*parm_.delta) )+ (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta))/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[i]-I*parm_.delta) );
                    
                    BCD[1] += (Ham_.VX_LR[i*M+j]*Ham_.VY_LR[j*M+i]-Ham_.VY_LR[i*M+j]*Ham_.VX_LR[j*M+i])*Ham_.VY_LR[i*(M+1)] * FD(w,parm_.T) *((1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta))*(1.0/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[j] + I*parm_.delta) + 1.0/(w-parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-parm_.W-Ham_.EN[j] + I*parm_.delta) )+ (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta))/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[i]-I*parm_.delta) );

                    BCD[2] += (Ham_.VZ_LR[i*M+j]*Ham_.VX_LR[j*M+i]-Ham_.VX_LR[i*M+j]*Ham_.VZ_LR[j*M+i])*Ham_.VZ_LR[i*(M+1)] * FD(w,parm_.T) *((1.0/(w-Ham_.EN[i]+I*parm_.delta)-1.0/(w-Ham_.EN[i]-I*parm_.delta))*(1.0/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[j] + I*parm_.delta) + 1.0/(w-parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-parm_.W-Ham_.EN[j] + I*parm_.delta) )+ (1.0/(w-Ham_.EN[j]+I*parm_.delta)-1.0/(w-Ham_.EN[j]-I*parm_.delta))/(w+parm_.W-Ham_.EN[i] + I*parm_.delta)/(w-Ham_.EN[i]-I*parm_.delta) );
                }
            }
        }

        //DrudeL_ += dw * real(DrudeL)/(2.0*pi);
        for (int i = 0; i < 3; i++){
            Inj_[i] += dw * imag(Inj[i])/(2.0*pi);
            Inj2_[i] += dw * imag(Inj2[i])/(2.0*pi);
            BCD_[i] += dw * imag(BCD[i])/(2.0*pi);
        }
    } 
};


