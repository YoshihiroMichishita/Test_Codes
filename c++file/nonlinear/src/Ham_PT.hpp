#pragma once

#include "const.hpp"
#include "parm_PT.hpp"
#include "matrix_op_mypc.hpp"


using namespace std;

// THis is the class for calculating the concrete Hamiltonian, velocity operator, and Green functions.

static const int M = 4;
static const int Mf = 4;
static const int D = 2;



class Ham{
    public:
        Complex H_k[M*M];
        Complex VX[M*M],VY[M*M],VZ[M*M];
        Complex VXX[M*M],VYX[M*M],VYY[M*M],VYXX[M*M],VYYX[M*M];
        Complex VR_k[M*M],VL_b[M*M],VR_b[M*M],VL_k[M*M],E_NH[M];
        double EN[M];
        Complex VX_LL[M*M],VX_LR[M*M],VX_RR[M*M],VZ_LL[M*M],VZ_LR[M*M],VZ_RR[M*M],VY_RR[M*M],VY_LR[M*M],VY_LL[M*M];
        Complex VXX_LR[M*M],VXX_LL[M*M],VXX_RR[M*M],VYX_LL[M*M],VYX_LR[M*M],VYX_RR[M*M]
                    ,VYY_LR[M*M],VYY_RR[M*M],VYY_LL[M*M],VYXX_LR[M*M],VYYX_LR[M*M];
        double NH_fac[M];
        Ham();
        Ham(parm parm_,double k[D],int sw);
        ~Ham();
        void Ham_list(int sw);
        void H_mom_NH(parm parm_, double k[D]);
        void BI_Velocity();
        void NH_factor();
};


class Green{
    public:
        Complex GR[M*M],dGR[M*M],GA[M*M],GRp[M*M],GRm[M*M],GAp[M*M],GAm[M*M],GRpp[M*M],GAmm[M*M],GRmA[M*M];
        Green();
        Green(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M]);
        Green(parm parm_,double w, double dw, double im[Mf],double re[Mf],Complex H[M*M]);
        ~Green();
        void G_List();
};
        
void H_mom(parm parm_, double k[D], Complex H[M*M]);
void H_mom_BI(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VL_b[M*M],Complex E[M]);
//void BI_Velocity(Ham Ham_);
//void H_mom_NH(parm parm_, double k[D],Ham Ham_);
//void H_mom_NH(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M],Complex E_NH[M]);
//void NH_factor(Complex Vx[M*M],Complex Vxx[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M]
//    , Complex Vx_LL[M*M], Complex Vx_RR[M*M], Complex Vx_LR[M*M],Complex Vxx_LL[M*M],Complex Vxx_LR[M*M], double NH_fac[M]);
//void NH_factor(Ham Ham_);
void Vx(parm parm_, double k[D], Complex H[M*M]);
void Vy(parm parm_, double k[D], Complex H[M*M]);
void Vz(parm parm_, double k[D], Complex H[M*M]);
void Vxx(parm parm_, double k[D], Complex H[M*M]);
void Vyx(parm parm_, double k[D], Complex H[M*M]);
void Vyxx(parm parm_, double k[D], Complex H[M*M]);
void Vyy(parm parm_, double k[D], Complex H[M*M]);
void Vyyx(parm parm_, double k[D], Complex H[M*M]);
void GreenR_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]);
void dGreenR_mom(parm parm_, double w,double dw, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M], Complex dG[M*M]);
void dGreenR_mom2(Complex G[M*M], Complex dG[M*M]);
void GreenA_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]);
void GreenR_minusA(Complex GR[M*M], Complex GA[M*M], Complex G[M*M]);
void GreenR_mom_p(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M],Complex GRp[M*M]);
void GreenR_mom_m(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M],Complex GRm[M*M]);
void GreenA_mom_p(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M],Complex GAp[M*M]);
void GreenA_mom_m(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M],Complex GAm[M*M]);
void GreenR_mom_pp(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GRpp[M*M]);
void GreenA_mom_mm(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GAmm[M*M]);

