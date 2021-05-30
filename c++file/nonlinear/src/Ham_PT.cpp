#include "Ham_PT.hpp"
#define s sin
#define c cos

using namespace std;


double a1[D];
double a2[D];

void H_mom(parm parm_, double k[D], Complex H[M*M]){

    double s_x = parm_.hx + (-parm_.a_R + parm_.a_D) * s(k[0]);
    double s_y = parm_.hy + (parm_.a_R + parm_.a_D) * s(k[1]);
    double s_z = parm_.hz;
    double s_0 = -parm_.t_i * (c(k[0]) + c(k[1]))+parm_.mu;
    double Vk = -parm_.t_e * c(k[0]/2) * c(k[1]/2);

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

}

void H_mom_BI(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VL_b[M*M], double E[M]){
    double s_x = parm_.hx + (-parm_.a_R + parm_.a_D) * s(k[0]);
    double s_y = parm_.hy + (parm_.a_R + parm_.a_D) * s(k[1]);
    double s_z = parm_.hz;
    double s_0 = -parm_.t_i * (c(k[0]) + c(k[1]))+parm_.mu;
    double Vk = -parm_.t_e * c(k[0]/2) * c(k[1]/2);

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

    Diag_H<M>(H,VR_k,E);
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
        VL_b[i*M+j] = conj(VR_k[j*M+i]);
        }
    }
};

/*
void BI_Velocity(Ham Ham_){
    Prod3<M>(Ham_.VL_b,Ham_.VX,Ham_.VR_k,Ham_.VX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VY,Ham_.VR_k,Ham_.VY_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VXX,Ham_.VR_k,Ham_.VXX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYX,Ham_.VR_k,Ham_.VYX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYY,Ham_.VR_k,Ham_.VYY_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYXX,Ham_.VR_k,Ham_.VYXX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYYX,Ham_.VR_k,Ham_.VYYX_LR);
};*/


/*
void NH_factor(Complex Vx[M*M],Complex Vxx[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M]
    , Complex Vx_LL[M*M], Complex Vx_RR[M*M], Complex Vx_LR[M*M],Complex Vxx_LL[M*M],Complex Vxx_LR[M*M], double NH_fac[M]){
    
    Prod3<M>(VL_b,Vx,VR_k,Vx_LR);
    Prod3<M>(VR_b,Vx,VR_k,Vx_RR);
    Prod3<M>(VL_b,Vx,VL_k,Vx_LL);
    Prod3<M>(VL_b,Vy,VR_k,Vy_LR);
    Prod3<M>(VR_b,Vy,VR_k,Vy_RR);
    Prod3<M>(VL_b,Vy,VL_k,Vy_LL);
    Prod3<M>(VL_b,Vxx,VL_k,Vxx_LL);
    Prod3<M>(VL_b,Vxx,VR_k,Vxx_LR);

    for (int i = 0; i < M; i++){
        NH_fac[i] = real((VR_b[i*M]*VR_k[i] + VR_b[i*M+1]*VR_k[i+M]) * (VL_b[i*M]*VL_k[i] + VL_b[i*M+1]*VL_k[i+M]));
    }
    
};
*/
/*
void NH_factor(Ham Ham_){
    
    Prod3<M>(Ham_.VL_b,Ham_.VX,Ham_.VR_k,Ham_.VX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VX,Ham_.VL_k,Ham_.VX_LL);
    Prod3<M>(Ham_.VR_b,Ham_.VX,Ham_.VR_k,Ham_.VX_RR);

    Prod3<M>(Ham_.VL_b,Ham_.VY,Ham_.VR_k,Ham_.VY_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VY,Ham_.VL_k,Ham_.VY_LL);
    Prod3<M>(Ham_.VR_b,Ham_.VY,Ham_.VR_k,Ham_.VY_RR);

    Prod3<M>(Ham_.VL_b,Ham_.VXX,Ham_.VR_k,Ham_.VXX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VXX,Ham_.VL_k,Ham_.VXX_LL);
    Prod3<M>(Ham_.VR_b,Ham_.VXX,Ham_.VR_k,Ham_.VXX_RR);

    Prod3<M>(Ham_.VL_b,Ham_.VYX,Ham_.VR_k,Ham_.VYX_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYX,Ham_.VL_k,Ham_.VYX_LL);
    Prod3<M>(Ham_.VR_b,Ham_.VYX,Ham_.VR_k,Ham_.VYX_RR);

    Prod3<M>(Ham_.VL_b,Ham_.VYY,Ham_.VR_k,Ham_.VYY_LR);
    Prod3<M>(Ham_.VL_b,Ham_.VYY,Ham_.VL_k,Ham_.VYY_LL);
    Prod3<M>(Ham_.VR_b,Ham_.VYY,Ham_.VR_k,Ham_.VYY_RR);

    for (int i = 0; i < M; i++){
        Ham_.NH_fac[i] = real((Ham_.VR_b[i*M]*Ham_.VR_k[i] + Ham_.VR_b[i*M+1]*Ham_.VR_k[i+M])
                                 * (Ham_.VL_b[i*M]*Ham_.VL_k[i] + Ham_.VL_b[i*M+1]*Ham_.VL_k[i+M]));
    }
    
};*/


void Vx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = (-parm_.a_R + parm_.a_D) * c(k[0]);
    double s_y = 0;
    double s_z = 0;
    double s_0 = parm_.t_i * s(k[0]);
    double Vk = parm_.t_e * s(k[0]/2) * c(k[1]/2)/2;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

}

void Vy(parm parm_, double k[D],Complex H[M*M]){
    double s_x = 0;
    double s_y = (parm_.a_R + parm_.a_D) * c(k[1]);
    double s_z = 0;
    double s_0 = parm_.t_i * s(k[1]);
    double Vk = parm_.t_e * c(k[0]/2) * s(k[1]/2)/2;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;
}

void Vz(parm parm_, double k[D], Complex H[M*M]){
    double s_x = 0;
    double s_y = 0;
    double s_z = 0;
    double s_0 = 0;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y;
    H[2] = s_x + I * s_y; H[3] = -s_z + s_0;
}

void Vyx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = (-parm_.a_R + parm_.a_D) * c(k[0]);
    double s_y = 0;
    double s_z = 0;
    double s_0 = 0;
    double Vk = -parm_.t_e * s(k[0]/2) * s(k[1]/2)/4;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

}

void Vxx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = -(-parm_.a_R + parm_.a_D) * s(k[0]);
    double s_y = 0;
    double s_z = 0;
    double s_0 = parm_.t_i * c(k[0]);
    double Vk = parm_.t_e * c(k[0]/2) * c(k[1]/2)/4;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

}

void Vyxx(parm parm_, double k[D], Complex H[M*M]){

    double s_x = 0;
    double s_y = 0;
    double s_z = 0;
    double s_0 = 0;
    double Vk = -parm_.t_e * c(k[0]/2) * s(k[1]/2)/8;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;

}

void Vyy(parm parm_, double k[D],Complex H[M*M]){
    double s_x = 0;
    double s_y = -(parm_.a_R + parm_.a_D) * s(k[1]);
    double s_z = 0;
    double s_0 = parm_.t_i * c(k[1]);
    double Vk = parm_.t_e * c(k[0]/2) * c(k[1]/2)/4;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;
}

void Vyyx(parm parm_, double k[D],Complex H[M*M]){
    double s_x = 0;
    double s_y = 0;
    double s_z = 0;
    double s_0 = 0;
    double Vk = -parm_.t_e * s(k[0]/2) * c(k[1]/2)/8;

    
    H[0] = s_z + s_0; H[1] = s_x - I * s_y; H[2] = Vk; H[3] = 0;
    H[4] = s_x + I * s_y; H[5] = -s_z + s_0; H[6] = 0; H[7] = Vk;
    H[8] = Vk; H[9] = 0; H[10] = (-s_z+s_0); H[11] = -(s_x - I * s_y);
    H[12] = 0; H[13] = Vk; H[14] = -(s_x + I * s_y); H[15] = s_z + s_0;
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
    G0[0] += -I * parm_.NH; G0[3] += I * parm_.NH;
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
    G0[0] += -I * parm_.NH; G0[3] += I * parm_.NH;
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
    G0[0] += I * parm_.NH; G0[3] += -I * parm_.NH;
    inverse_NH<M>(G0,G);
    
}

void GreenR_minusA(Complex GR[M*M], Complex GA[M*M], Complex G[M*M]){
    for (int i = 0; i < M*M; i++){
        G[i] = GR[i]-GA[i];
    }
};

void GreenR_mom_p(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M],Complex GRp[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w+parm_.W) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    G0[0] += -I * parm_.NH; G0[3] += I * parm_.NH;
    inverse_NH<M>(G0,GRp);
};

void GreenR_mom_m(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M], Complex GRm[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w-parm_.W) + I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] - I * im[Mf-i-1]; 
    }
    G0[0] += -I * parm_.NH; G0[3] += I * parm_.NH;
    inverse_NH<M>(G0,GRm);
};

void GreenA_mom_p(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M], Complex GAp[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w+parm_.W) - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    G0[0] += I * parm_.NH; G0[3] += -I * parm_.NH;
    inverse_NH<M>(G0,GAp);
};

void GreenA_mom_m(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M], Complex GAm[M*M]){
    Complex G0[M*M];
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            G0[i*M+j] = -H[i*M+j];
        }
        G0[i*(M+1)] += parm_.alpha * (w-parm_.W) - I * parm_.delta;    
    }

    for (int i = 0; i < Mf; i++){
        G0[(M-i)*M-1-i] += - re[Mf-i-1] + I * im[Mf-i-1]; 
    }
    G0[0] += I * parm_.NH; G0[3] += -I * parm_.NH;
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
    G0[0] += -I * parm_.NH; G0[3] += I * parm_.NH;
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
    G0[0] += I * parm_.NH; G0[3] += -I * parm_.NH;
    inverse_NH<M>(G0,GAmm);
};

Ham::Ham(){
};

Ham::Ham(parm parm_,double k[D],int sw){
    
    //calculate with the conventional band index
    if(sw==0){
        H_mom_BI(parm_,k,H_k,VR_k,VL_b,EN);
        Vx(parm_,k,VX);
        Vy(parm_,k,VY);
        //Vz(parm_,k,VZ);
        Vxx(parm_,k,VXX);
        Vyxx(parm_,k,VYXX);
        Vyx(parm_,k,VYX);
        Vyy(parm_,k,VYY);
        Vyyx(parm_,k,VYYX);
        BI_Velocity();
    }

    //calculate with the Green function method
    else if(sw==1){
        H_mom(parm_,k,H_k);
        Vx(parm_,k,VX);
        Vy(parm_,k,VY);
        //Vz(parm_,k,VZ);
        Vxx(parm_,k,VXX);
        Vy(parm_,k,VY);
        Vyx(parm_,k,VYX);
        //Vyxx(parm_,k,VYXX);
    }

    //calculate with the Non-Hermitian band index
    else if(sw==2){
        //H_mom_NH(parm_,k,H_k,VR_k,VR_b,VL_k,VL_b,E_NH);
        H_mom_NH(parm_,k);
        Vx(parm_,k,VX);
        Vxx(parm_,k,VXX);
        Vy(parm_,k,VY);
        Vyx(parm_,k,VYX);
        Vyxx(parm_,k,VYXX);
        NH_factor();
    }
    else{
        cout << "error!" << endl; 
    }
};


Ham::~Ham(){
};

void Ham::Ham_list(int sw){
    cout << "=================" << endl;
    cout << "Hamiltonian" << endl;
    cout << "--------------" << endl;
    for (int i = 0; i < M; i++){
        cout << "( " ;
        for (int j = 0; j < M; j++){
            cout << H_k[i*M+j] << ", ";
        }
        cout << ") " << endl;
    }
    if(sw>0){
        for (int i = 0; i < M; i++){
            cout << "*****************" << endl;
            cout << "EigenVector " << i << " " << EN[i] <<  endl;
            cout << "--------------" << endl;
            cout << "( " ;
            Complex test[M];
            for (int j = 0; j < M; j++){
                cout << VR_k[j*M+i] << ", ";
                for (int l = 0; l < M; l++){
                    test[j] += H_k[j*M+l]*VR_k[l*M+i];
                }
                
            }
            cout << ")  H|i> = (";
            for (int j = 0; j < M; j++){
                cout << test[j] << ", ";
            }
            cout << ")" << endl;
        }
        if(sw == 2){
            for (int i = 0; i < M; i++){
                cout << "*****************" << endl;
                cout << "EigenVector(b) " << i << " " << EN[i] <<  endl;
                cout << "--------------" << endl;
                cout << "( " ;
                Complex test[M];
                for (int j = 0; j < M; j++){
                    cout << VL_b[i*M+j] << ", ";
                    for (int l = 0; l < M; l++){
                        test[j] += VL_b[i*M+l]*H_k[l*M+j];
                    }
                }
                cout << ")  <i|H = (";
                for (int j = 0; j < M; j++){
                    cout << test[j] << ", ";
                }
                cout << ")" << endl;
            }
        }
        cout << "*****************" << endl;
        cout << "Vx " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VX[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vx_BI " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VX_LR[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vy " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VY[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vy_BI " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VY_LR[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vxx " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VXX[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vxx_BI " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VXX_LR[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vyx " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VYX[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
        cout << "*****************" << endl;
        cout << "Vyx_BI " << endl;
        cout << "--------------" << endl;
        for (int i = 0; i < M; i++){
            cout << "( " ;
            for (int j = 0; j < M; j++){
                cout << VYX_LR[i*M+j] << ", ";
            }
            cout << ") " << endl;
        }
    }
    
}


void Ham::H_mom_NH(parm parm_, double k[D]){

    double s_x = parm_.hx + (-parm_.a_R + parm_.a_D) * s(k[0]);
    double s_y = parm_.hy + (parm_.a_R + parm_.a_D) * s(k[1]);
    double s_z = parm_.hz;
    double s_0 = -parm_.t_i * (c(k[0]) + c(k[1]));
    double Vk = -parm_.t_e * c(k[0]/2) * c(k[1]/2);

    
    H_k[0] = s_z + s_0 - I*parm_.delta; H_k[1] = s_x - I * s_y; H_k[2] = Vk; H_k[3] = 0;
    H_k[4] = s_x + I * s_y; H_k[5] = -s_z + s_0 - I*parm_.delta; H_k[6] = 0; H_k[7] = Vk;
    H_k[8] = Vk; H_k[9] = 0; H_k[10] = (-s_z+s_0) - I*parm_.delta; H_k[11] = -(s_x - I * s_y);
    H_k[12] = 0; H_k[13] = Vk; H_k[14] = -(s_x + I * s_y); H_k[15] = s_z + s_0 - I*parm_.delta;

    Complex G[M*M];
    for (int i = 0; i < M*M; i++){
        G[i] = H_k[i];
    }
    

    //Diag_NH<M>(G,VL_k,VR_k,E_NH);
    Diag_NH<M>(G,VL_b,VR_k,E_NH);

    
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


void Ham::NH_factor(){
    
    Prod3<M>(VL_b,VX,VR_k,VX_LR);
    Prod3<M>(VL_b,VX,VL_k,VX_LL);
    Prod3<M>(VR_b,VX,VR_k,VX_RR);

    Prod3<M>(VL_b,VY,VR_k,VY_LR);
    Prod3<M>(VL_b,VY,VL_k,VY_LL);
    Prod3<M>(VR_b,VY,VR_k,VY_RR);

    Prod3<M>(VL_b,VXX,VR_k,VXX_LR);
    Prod3<M>(VL_b,VXX,VL_k,VXX_LL);
    Prod3<M>(VR_b,VXX,VR_k,VXX_RR);

    Prod3<M>(VL_b,VYX,VR_k,VYX_LR);
    Prod3<M>(VL_b,VYX,VL_k,VYX_LL);
    Prod3<M>(VR_b,VYX,VR_k,VYX_RR);

    Prod3<M>(VL_b,VYY,VR_k,VYY_LR);
    Prod3<M>(VL_b,VYY,VL_k,VYY_LL);
    Prod3<M>(VR_b, VYY, VR_k, VYY_RR);

    for (int i = 0; i < M; i++){
        NH_fac[i] = real((VR_b[i*M]*VR_k[i] + VR_b[i*M+1]*VR_k[i+M])
                                 * (VL_b[i*M]*VL_k[i] + VL_b[i*M+1]*VL_k[i+M]));
    }
    
};


void Ham::BI_Velocity(){
    Prod3<M>(VL_b,VX,VR_k,VX_LR);
    Prod3<M>(VL_b,VY,VR_k,VY_LR);
    Prod3<M>(VL_b,VXX,VR_k,VXX_LR);
    Prod3<M>(VL_b,VYX,VR_k,VYX_LR);
    Prod3<M>(VL_b,VYY,VR_k,VYY_LR);
    Prod3<M>(VL_b,VYXX,VR_k,VYXX_LR);
    Prod3<M>(VL_b,VYYX,VR_k,VYYX_LR);
};


Green::Green(){
};

Green::Green(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M]){
    GreenR_mom(parm_,w,im,re,H,GR);
    GreenA_mom(parm_,w,im,re,H,GA);
    GreenR_minusA(GR,GA,GRmA);
    dGreenR_mom2(parm_,GR,dGR);
    GreenR_mom_p(parm_,w,im,re,H,GRp);
    GreenR_mom_m(parm_,w,im,re,H,GRm);
    GreenA_mom_p(parm_,w,im,re,H,GAp);
    GreenA_mom_m(parm_,w,im,re,H,GAm);
};

Green::Green(parm parm_,double w, double dw, double im[Mf],double re[Mf],Complex H[M*M]){
    GreenR_mom(parm_,w,im,re,H,GR);
    GreenA_mom(parm_,w,im,re,H,GA);
    GreenR_minusA(GR,GA,GRmA);
    //dGreenR_mom(parm_,w,dw,im,re,H,GR,dGR);
    dGreenR_mom2(parm_,GR,dGR);
    GreenR_mom_p(parm_,w,im,re,H,GRp);
    GreenR_mom_m(parm_,w,im,re,H,GRm);
    GreenA_mom_p(parm_,w,im,re,H,GAp);
    GreenA_mom_m(parm_,w,im,re,H,GAm);
    GreenR_mom_pp(parm_,w,parm_.W,im,re,H,GRpp);
    GreenA_mom_mm(parm_,w,parm_.W,im,re,H,GAmm);
};


Green::~Green(){
};

void Green::G_List(){

    cout << "=============" << endl;
    cout << "Green R: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GR[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green A: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GA[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green R-A: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GRmA[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green Rp: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GRp[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green Rm: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GRm[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green Ap: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GAp[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    cout << "=============" << endl;
    cout << "Green Am: " << endl;
    for (int i = 0; i < M; i++){
        cout << "(";
        for (int j = 0; j < M; j++){
            cout << GAm[i*M+j] << ", ";
        }
        cout << ")" << endl;
    }
    
}