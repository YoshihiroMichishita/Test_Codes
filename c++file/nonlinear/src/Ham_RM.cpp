#include "Ham_RM.hpp"
#define s sin
#define c cos

using namespace std;


double a1[D];
double a2[D];

void H_mom(parm parm_, double k[D], Complex H[M*M]){

    double s_x = parm_.t_i * c(k[0]);
    double s_y = parm_.t_e * s(k[0]);
    double s_z = parm_.a_d;

    
    H[0] = s_z + parm_.mu; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z + parm_.mu;

}

void H_mom_BI(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VL_b[M*M], double E[M]){
    double s_x = parm_.t_i * c(k[0]);
    double s_y = parm_.t_e * s(k[0]);
    double s_z = parm_.a_d;

    
    H[0] = s_z + parm_.mu; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z + parm_.mu;

    Diag_H<M>(H,VR_k,E);
    for (int i = 0; i < M*M; i++){
        VL_b[i] = conj(VR_k[i]);
    }
};

void BI_Velocity(Complex Vx[M*M],Complex Vxx[M*M], Complex VR_k[M*M], Complex VL_b[M*M], Complex Vx_LR[M*M], Complex Vxx_LR[M*M]){
    Prod3<M>(VL_b,Vx,VR_k,Vx_LR);
    //InPro_M<M>(VL_b,Vy,VR_k,Vy_LR);
    Prod3<M>(VL_b,Vxx,VR_k,Vxx_LR);
    //InPro_M<M>(VL_b,Vyx,VR_k,Vyx_LR);
};

void H_mom_NH(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M],Complex E_NH[M]){

    double s_x = parm_.t_i * c(k[0]);
    double s_y = parm_.t_e * s(k[0]);
    Complex s_z = parm_.a_d + I*parm_.NH;

    
    H[0] = s_z + parm_.mu; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z + parm_.mu;

    Complex G[M*M];
    G[0] = H[0]; G[1] = H[1]; 
    G[2] = H[2]; G[3] = H[3];

    //Diag_NH<M>(G,VL_k,VR_k,E_NH);
    Diag_NH<M>(G,VL_b,VR_k,E_NH);

    /*
    Complex VL_k2[M*M];
    for (int i = 0; i < M; i++){
        VL_k2[i] = VL_k[i]/(conj(VL_k[i])*VR_k[i] + conj(VL_k[i+M])*VR_k[i+M]);
        VL_k2[i+M] = VL_k[i+M]/(conj(VL_k[i])*VR_k[i] + conj(VL_k[i+M])*VR_k[i+M]);
    }*/

    
    Complex VL_b2[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            VL_b2[i*M+j] = VL_b[i*M+j]/(VL_b[i*M]*VR_k[i] + VL_b[i*M+1]*VR_k[i+M]);
        }
    }
    for (int i = 0; i < M*M; i++){
        VL_b[i] = VL_b2[i];
    }
    
    
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            VR_b[i*M+j] = conj(VR_k[i+M*j]);
            VL_k[i*M+j] = conj(VL_b[i+M*j]);
        }
    }
}

void NH_factor(Complex Vx[M*M],Complex Vxx[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M]
    , Complex Vx_LL[M*M], Complex Vx_RR[M*M], Complex Vx_LR[M*M],Complex Vxx_LL[M*M],Complex Vxx_LR[M*M], double NH_fac[M]){
    
    /*
    InPro_M<M>(VL_b,Vx,VR_k,Vx_LR);
    InPro_M<M>(VR_b,Vx,VR_k,Vx_RR);
    InPro_M<M>(VL_b,Vx,VL_k,Vx_LL);
    InPro_M<M>(VR_b,Vy,VR_k,Vy_RR);
    InPro_M<M>(VL_b,Vy,VR_k,Vy_LR);
    InPro_M<M>(VL_b,Vxx,VL_k,Vxx_LL);
    InPro_M<M>(VL_b,Vxx,VR_k,Vxx_LR);*/
    Prod3<M>(VL_b,Vx,VR_k,Vx_LR);
    Prod3<M>(VR_b,Vx,VR_k,Vx_RR);
    Prod3<M>(VL_b,Vx,VL_k,Vx_LL);
    Prod3<M>(VL_b,Vxx,VL_k,Vxx_LL);
    Prod3<M>(VL_b,Vxx,VR_k,Vxx_LR);

    for (int i = 0; i < M; i++){
        NH_fac[i] = real((VR_b[i*M]*VR_k[i] + VR_b[i*M+1]*VR_k[i+M]) * (VL_b[i*M]*VL_k[i] + VL_b[i*M+1]*VL_k[i+M]));
    }
    
};

void Vx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = -parm_.t_i * s(k[0]);
    double s_y = parm_.t_e * c(k[0]);
    double s_z = 0;

    
    H[0] = s_z; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z;   

}

void Vxx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = -parm_.t_i * c(k[0]);
    double s_y = -parm_.t_e * s(k[0]);
    double s_z = 0;

    
    H[0] = s_z; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z; 

}

void GreenR_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * w + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,G);
}

void dGreenR_mom(parm parm_, double w,double dw, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M], Complex dG[M*M]){
    Complex G0[M*M],G1[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w + dw) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,G1);
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            dG[i*M+j] = (G1[i*M+j]-G[i*M+j])/dw;
        }   
    }
}

void dGreenR_mom2(parm parm_,Complex G[M*M], Complex dG[M*M]){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            for (int k = 0; k < M; k++){
                dG[i*M+j] += -parm_.alpha*G[i*M+k]*G[k*M+j];
            }
        }   
    }
};

void GreenA_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * w - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,G);
    
}

void GreenR_minusA(Complex GR[M*M], Complex GA[M*M], Complex G[M*M]){
    for (int i = 0; i < M*M; i++){
        G[i] = GR[i]-GA[i];
    }
};

void GreenR_mom_p(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GRp[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w+W) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GRp);
};

void GreenR_mom_m(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M], Complex GRm[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w-W) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GRm);
};

void GreenA_mom_p(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M], Complex GAp[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w+W) - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GAp);
};

void GreenA_mom_m(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M], Complex GAm[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w-W) - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GAm);
};

void GreenR_mom_pp(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M], Complex GRpp[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w+2*W) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GRpp);
};

void GreenA_mom_mm(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M], Complex GAmm[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w-2*W) - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    inverse_NH<M>(G0,GAmm);
};

Ham::Ham(){
};

Ham::Ham(parm parm_,double k[D],int sw){
    
    //calculate with the conventional band index
    if(sw==0){
        H_mom_BI(parm_,k,H_k,VR_k,VL_b,EN);
        Vx(parm_,k,VX);
        Vxx(parm_,k,VXX);
        //Vy(parm_,k,VY);
        //Vyx(parm_,k,VYX);
        //Vyxx(parm_,k,VYXX);
        BI_Velocity(VX,VXX,VR_k,VL_b,VX_LR,VXX_LR);
    }

    //calculate with the Green function method
    else if(sw==1){
        H_mom(parm_,k,H_k);
        Vx(parm_,k,VX);
        Vxx(parm_,k,VXX);
        //Vy(parm_,k,VY);
        //Vyx(parm_,k,VYX);
        //Vyxx(parm_,k,VYXX);
    }

    //calculate with the Non-Hermitian band index
    else if(sw==2){
        H_mom_NH(parm_,k,H_k,VR_k,VR_b,VL_k,VL_b,E_NH);
        Vx(parm_,k,VX);
        Vxx(parm_,k,VXX);
        //Vy(parm_,k,VY);
        //Vyx(parm_,k,VYX);
        //Vyxx(parm_,k,VYXX);
        NH_factor(VX,VXX,VR_k,VR_b,VL_k,VL_b,VX_LL,VX_RR,VX_LR,VXX_LL,VXX_LR,NH_fac);
    }
    else{
        cout << "error!" << endl; 
    }
};


Ham::~Ham(){
};

Green::Green(){
};

Green::Green(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M]){
    GreenR_mom(parm_,w,im,re,H,GR);
    GreenA_mom(parm_,w,im,re,H,GA);
    dGreenR_mom2(parm_,GR,dGR);
};

Green::Green(parm parm_,double w, double dw, double im[Mf],double re[Mf],Complex H[M*M]){
    GreenR_mom(parm_,w,im,re,H,GR);
    GreenA_mom(parm_,w,im,re,H,GA);
    GreenR_minusA(GR,GA,GRmA);
    //dGreenR_mom(parm_,w,dw,im,re,H,GR,dGR);
    dGreenR_mom2(parm_,GR,dGR);
    GreenR_mom_p(parm_,w,parm_.W,im,re,H,GRp);
    GreenR_mom_m(parm_,w,parm_.W,im,re,H,GRm);
    GreenA_mom_p(parm_,w,parm_.W,im,re,H,GAp);
    GreenA_mom_m(parm_,w,parm_.W,im,re,H,GAm);
    GreenR_mom_pp(parm_,w,parm_.W,im,re,H,GRpp);
    GreenA_mom_mm(parm_,w,parm_.W,im,re,H,GAmm);
};


Green::~Green(){
};
