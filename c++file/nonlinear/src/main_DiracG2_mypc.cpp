#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out1("Disp.dat");
    ofstream out2("Cond_kz.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    double dw = 2.0*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = 4.0 * pi *pi *pi / (parm_.K_SIZE*parm_.K_SIZE*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;

    double DrudeL_sum[parm_.W_SIZE],Inj_sum[parm_.W_SIZE][3],Inj2_sum[parm_.W_SIZE][3];
    double DrudeL_total,BCD_total[3],Inj_total[3],Inj2_total[3];
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    #pragma omp parallel for num_threads(52)
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        double w = 2.0*(ww-parm_.W_SIZE/2) * parm_.W_MAX /parm_.W_SIZE;
        double DrudeL_=0,Inj_[3]={0,0,0},Inj2_[3]={0,0,0};

        for (int i = 0; i < parm_.K_SIZE; i++){
            double kz = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
            for (int j = 0; j < parm_.K_SIZE; j++){
                double kx = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
                for (int l = 0; l < parm_.K_SIZE; l++){
                    double ky = 2.0*(l-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
                    if(fabs(kz)+fabs(kx)+fabs(ky)<=(3*pi/2)){
                        double k[3] = {kx,ky,kz};
                        Ham Ham_(parm_,k,0);
                        Opt_Green_transport_BI2(parm_,Ham_,w,dw,DrudeL_,Inj_,Inj2_);
                    }
                }
            }
            
        }
        DrudeL_sum[ww] = dk2 * DrudeL_ / (8*pi*pi*pi);
        for (int j = 0; j < 3; j++){
            Inj_sum[ww][j] = dk2 * Inj_[j] / (8*pi*pi*pi);
            Inj2_sum[ww][j] = dk2 * Inj2_[j] / (8*pi*pi*pi);
        }
    }
    
        

        
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        double w = 2.0*(i-parm_.W_SIZE/2)*parm_.W_MAX/parm_.W_SIZE;
        DrudeL_total += DrudeL_sum[i];
        out2 << w;

        for (int j = 0; j < 3; j++){
            if(i==0 || i==parm_.W_SIZE){
                Inj_total[j] += Inj_sum[i][j]/3;
                Inj2_total[j] += Inj2_sum[i][j]/3;
            }
            else if(i%2==1){
                Inj_total[j] += 4*Inj_sum[i][j]/3;
                Inj2_total[j] += 4*Inj2_sum[i][j]/3;
            }
            else{
                Inj_total[j] += 2*Inj_sum[i][j]/3;
                Inj2_total[j] += 2*Inj2_sum[i][j]/3;
            }
            out2 << "\t" << Inj_sum[i][j] << "\t" << Inj2_sum[i][j];
        }

        out2 << endl;
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    cout << "DrudeL " << DrudeL_total << endl;
    cout << "Inj " << Inj_total[0] << " " << Inj_total[1] << " " << Inj_total[2] << endl;
    cout << "Inj2 " << Inj2_total[0] << " " << Inj2_total[1] << " " << Inj2_total[2] << endl;
    //cout << "total " << Drude_sum+BCD_sum+Inj_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << parm_.delta << "\t" << DrudeL_total;
    for (int i = 0; i < 3; i++){
        out3  << "\t" << BCD_total[i] << "\t" << Inj_total[i] << "\t" << Inj2_total[i];
    }
    
    out3 << "\t" << Inj_total[0]+Inj_total[1]+Inj_total[2] << "\t" << Inj2_total[0]+Inj2_total[1]+Inj2_total[2] << endl;
    return 0;
}