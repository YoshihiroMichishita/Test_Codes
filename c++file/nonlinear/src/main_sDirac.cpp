#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out1("Disp.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    double dw = 2.0*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = 8.0 *parm_.K_MAX *parm_.K_MAX *parm_.K_MAX / (parm_.K_SIZE*parm_.K_SIZE*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;

    double DrudeL_sum[parm_.K_SIZE],BCD_sum[parm_.K_SIZE][3],Inj_sum[parm_.K_SIZE][3],Inj2_sum[parm_.K_SIZE][3];
   double DrudeL_total,BCD_total[3],Inj_total[3];
   //double Inj_total1=0,Inj_total2=0,Inj_total3=0;
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    #pragma omp parallel for// reduction(+:Inj_total1,Inj_total2,Inj_total3)
    for (int i = 0; i < parm_.K_SIZE; i++){

        double DrudeL_=0,Inj2_[3]={0,0,0},BCD_[3]={0,0,0},Inj_[3]={0,0,0};
        double kz = 2.0*(i-parm_.K_SIZE/2)*parm_.K_MAX/parm_.K_SIZE;

        for (int j = 0; j < parm_.K_SIZE; j++){
            double kx = 2.0*(j-parm_.K_SIZE/2)*parm_.K_MAX/parm_.K_SIZE;
            for (int l = 0; l < parm_.K_SIZE; l++){
                double ky = 2.0*(l-parm_.K_SIZE/2)*parm_.K_MAX/parm_.K_SIZE;
                //if(fabs(kz)+fabs(kx)+fabs(ky)<=(3*pi/2)){
                    double k[3] = {kx,ky,kz};
                    Ham Ham_(parm_,k,0);
                    //Opt_Green_transport_BI(parm_,Ham_,DrudeL_,BCD_,Inj_,Inj2_);
                    //Opt_RTA_transport_BI3(parm_,Ham_,DrudeL_,BCD_,Inj_);
                    Opt_RTA_transport_BI3(parm_,Ham_,Inj_);
                //}
                
            }
        }

        //Inj_total1 += dk2 * Inj_[0] / (8*pi*pi*pi);
        //Inj_total2 += dk2 * Inj_[1] / (8*pi*pi*pi);
        //Inj_total3 += dk2 * Inj_[2] / (8*pi*pi*pi);
        DrudeL_sum[i] = dk2 * DrudeL_ / (8*pi*pi*pi);
        for (int j = 0; j < 3; j++){
            BCD_sum[i][j] = dk2 * BCD_[j] / (8*pi*pi*pi);
            Inj_sum[i][j] = dk2 * Inj_[j] / (8*pi*pi*pi);
            //Inj2_sum[i][j] = dk2 * Inj2_[j] / (8*pi*pi*pi);
        }
    }
    
    
    
    for (int i = 0; i < parm_.K_SIZE; i++){
        DrudeL_total += DrudeL_sum[i];
        for (int j = 0; j < 3; j++){
            BCD_total[j] += BCD_sum[i][j];
            Inj_total[j] += Inj_sum[i][j];
            //Inj2_total[j] += Inj2_sum[i][j];
        }
        double kz = 2.0*(i-parm_.K_SIZE/2)*parm_.K_MAX/parm_.K_SIZE;
        double k2[3] = {0,0,kz};
        Ham Ham_(parm_,k2,0);
        out1 << kz << "\t" << Ham_.EN[0] << "\t" << Ham_.EN[1] << endl;
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    //cout << "DrudeL " << DrudeL_total << endl;
    //cout << "BCD " << BCD_total[0] << " " << BCD_total[1] << " " << BCD_total[2] << endl;
    cout << "Inj " << Inj_total[0] << " " << Inj_total[1] << " " << Inj_total[2] << endl;
    //cout << "Inj " << Inj_total1 << " " << Inj_total2 << " " << Inj_total3 << endl;
    //cout << "Inj2 " << Inj2_total[0] << " " << Inj2_total[1] << " " << Inj2_total[2] << endl;
    //cout << "total " << Drude_sum+BCD_sum+Inj_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << parm_.delta; //<< "\t" << Inj_total1 << "\t" << Inj_total2 << "\t" << Inj_total3 << "\t" << Inj_total1+Inj_total2+Inj_total3 << endl;
    
    for (int i = 0; i < 3; i++){
        out3  << "\t" << BCD_total[i] << "\t" << Inj_total[i];
    }
    
    out3 << "\t" << Inj_total[0]+Inj_total[1]+Inj_total[2] << endl;
    return 0;
}

