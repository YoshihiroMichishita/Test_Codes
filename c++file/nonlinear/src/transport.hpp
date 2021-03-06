#pragma once

//#include "Ham_TiltedDirac.hpp"
#include "Ham_PT.hpp"
//#include "Ham_TMD.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double FD(double w, double T);
double dFD(double w, double T);
double ddFD(double w, double T);
double Spectral(parm parm_,double w, double k[2], double re[Mf], double im[Mf]);

void DC_transport(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2);

void DC_transport_NRC(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& XXXr, double& XXXr2, double& XXXi, double& XXXi2);

void DC_NLH_NRC(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& YX, double& XXXr, double& XXXr2, double& XXXi, double& XXXi2, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2);

void Linear_transport_NH(double dw, double w,double delta,double T,Ham Ham_, double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi,
 double& YXXi2, double& XX_, double& YX_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_);

void Linear_transport_NHwithNRC(parm parm_,double dw, double w,double T,Ham Ham_, double& XX, double& XXX, double& XXX2, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2, double& XX_, double& XXX_, double& XXX2_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_);

void Linear_transport_NH_NRC(double dw, double w,double T,Ham Ham_, double& XX, double& XXX, double& XX_, double& XXX_);

void Opt_transport(parm parm_, double dw, double w,Ham Ham_,Green Green_,double& XX,double& XX2, double& PVYXi, double& div1, double& div2, double& div3);
void Opt_RTA_transport_BI(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_);
void Opt_RTA_transport_BI2(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_);
//void Opt_RTA_transport_BI3(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3]);
void Opt_RTA_transport_BI3(parm parm_,Ham Ham_,double Inj_[3]);
void Opt_RTA_transport_BI4(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3]);
void Opt_Green_transport_BI(parm parm_,Ham Ham_,double& DrudeL_, double BCD_[3], double Inj_[3] ,double Inj2_[3]);
void Opt_Green_transport_BI_div(parm parm_,Ham Ham_,double w,double dk2, double div_[3]);
void Opt_IFS_BI(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy);
void Opt_IFS_Green(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy);
void Opt_IFS_Green_w(parm parm_, Ham Ham_,double w, double& IFS_xxy, double& IFS_yxy);
void Opt_IFS_Green_w_div(parm parm_, Ham Ham_,double w, double& div_xxy, double& div_yxy);
void Opt_IFS_Green_test(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy);
void Opt_IFS_Green_test2(parm parm_, Ham Ham_, double& IFS_xxy, double& IFS_yxy);
void Opt_Gyration_BI(parm parm_, Ham Ham_, double& gyro_xxy, double& gyro_yxy);
void Opt_Gyration_BI2(parm parm_, Ham Ham_, double& gyro_xxy);
void Opt_Gyration_BI3(parm parm_, Ham Ham_, double& gyro_xxy, double& gyro_xxy_im);
void Opt_Gyration_Green(parm parm_, Ham Ham_, double w, double& gyro_xxy, double& gyro_yxy);
void Opt_Circular_Green(parm parm_, Ham Ham_,Green Green_,double w, double& PVCP1, double& PVCP2, double& PVCPWW, double& PVCP);
void Opt_Circular_Green_BI(parm parm_,Ham Ham_,double w, double& PVCP1, double& PVCP);