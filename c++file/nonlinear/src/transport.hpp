#pragma once

#include "Ham_TMD.hpp"
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

double FD(double w, double T);
double dFD(double w, double T);
double ddFD(double w, double T);
double Spectral(parm parm_,double w, double k[2], double re[Mf], double im[Mf]);

void DC_transport(double dw, double w,double T,Ham Ham_,Green Green_,double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi, double& YXXi2);

void Linear_transport_NH(double dw, double w,double T,Ham Ham_, double& XX, double& YX, double& YXXr, double& YXXr2, double& YXXi,
 double& YXXi2, double& XX_, double& YX_, double& YXXr_, double& YXXr2_, double& YXXi_, double& YXXi2_);

void Linear_transport_NH_NRC(double dw, double w,double T,Ham Ham_, double& XX, double& XXX, double& XX_, double& XXX_);

void Opt_transport(parm parm_, double dw, double w,Ham Ham_,Green Green_,double& XX,double& XX2, double& PVYXi, double& div1, double& div2, double& div3);
void Opt_RTA_transport_BI(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_);
void Opt_RTA_transport_BI2(parm parm_,Ham Ham_,double& DrudeL_ ,double& Drude_, double& BCD_, double& Inj_);
