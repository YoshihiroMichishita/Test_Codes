#pragma once

#include "const.hpp"
#include "parm_TMD.hpp"
#include "matrix_op_mypc.hpp"


using namespace std;

// This is the class for calculating the concrete Hamiltonian, velocity operator, and Green functions.

static const int M = 2;
static const int Mf = 2;
static const int D = 2;

        
void H_mom(parm parm_, double k[D], Complex H[M*M]);
void H_mom_BI(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VL_b[M*M],Complex E[M]);
void BI_Velocity(Complex Vx[M*M],Complex Vy[M*M],Complex Vxx[M*M],Complex Vyx[M*M], Complex VR_k[M*M], Complex VL_b[M*M], Complex Vx_LR[M*M], Complex Vy_LR[M*M], Complex Vxx_LR[M*M], Complex Vyx_LR[M*M]);
void H_mom_NH(parm parm_, double k[D], Complex H[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M],Complex E_NH[M]);
void NH_factor(Complex Vx[M*M],Complex Vy[M*M],Complex Vxx[M*M], Complex VR_k[M*M], Complex VR_b[M*M], Complex VL_k[M*M],Complex VL_b[M*M]
    , Complex Vx_LL[M*M], Complex Vx_RR[M*M], Complex Vx_LR[M*M], Complex Vy_RR[M*M], Complex Vy_LR[M*M],Complex Vxx_LL[M*M],Complex Vxx_LR[M*M], double NH_fac[M]);
void Vx(parm parm_, double k[D], Complex H[M*M]);
void Vxx(parm parm_, double k[D], Complex H[M*M]);
void Vy(parm parm_, double k[D], Complex H[M*M]);
void Vyx(parm parm_, double k[D], Complex H[M*M]);
void Vyxx(parm parm_, double k[D], Complex H[M*M]);
void GreenR_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]);
void dGreenR_mom(parm parm_, double w,double dw, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M], Complex dG[M*M]);
void dGreenR_mom2(parm parm_, Complex G[M*M], Complex dG[M*M]);
void GreenA_mom(parm parm_, double w, double im[Mf], double re[Mf],Complex H[M*M], Complex G[M*M]);
void GreenR_minusA(Complex GR[M*M], Complex GA[M*M], Complex G[M*M]);
void GreenR_mom_p(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GRp[M*M]);
void GreenR_mom_m(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GRm[M*M]);
void GreenA_mom_p(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GAp[M*M]);
void GreenA_mom_m(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GAm[M*M]);
void GreenR_mom_pp(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GRpp[M*M]);
void GreenA_mom_mm(parm parm_,double w,double W,double im[Mf],double re[Mf],Complex H[M*M],Complex GAmm[M*M]);


class Ham{
    public:
        Complex H_k[M*M];
        Complex VX[M*M],VY[M*M],VXX[M*M],VYX[M*M],VYXX[M*M];
        Complex VR_k[M*M],VL_b[M*M],VR_b[M*M],VL_k[M*M],E_NH[M];
        double EN[M];
        Complex VX_LL[M*M],VX_LR[M*M],VX_RR[M*M],VXX_LL[M*M],VXX_LR[M*M],VY_RR[M*M],VY_LR[M*M],VYX_LR[M*M],VZ_LR[M*M];
        double NH_fac[M];
        Ham();
        Ham(parm parm_, double k[D], int sw);
        ~Ham();
        void Ham_Check(int sw);
};

class Green{
    public:
        Complex GR[M*M],dGR[M*M],GA[M*M],GRp[M*M],GRm[M*M],GAp[M*M],GAm[M*M],GRpp[M*M],GAmm[M*M],GRmA[M*M];
        Green();
        Green(parm parm_,double w,double im[Mf],double re[Mf],Complex H[M*M]);
        Green(parm parm_,double w, double dw, double im[Mf],double re[Mf],Complex H[M*M]);
        ~Green();
};