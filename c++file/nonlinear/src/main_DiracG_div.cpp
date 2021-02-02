//#include "C3v_HSL.hpp"
#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out1("Div_w.dat");
    //ofstream out2("Cond_kz.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    double dw = 2.0*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = 8.0 * pi *pi *pi / (parm_.K_SIZE*parm_.K_SIZE*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;

    
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    int SIZE=0;
    for (int WW = 0; WW < parm_.W_SIZE; WW++){
        double w = 2.0*(WW-parm_.W_SIZE/2)*parm_.W_MAX/parm_.W_SIZE;
        if(w > 4*parm_.T){
            SIZE=WW;
            break;
        }
    }
    cout << "W_SIZE " << SIZE << endl;

    double div_sum[SIZE][3];
    double div_total[3];
    

    #pragma omp parallel for
    for (int WW = 0; WW < SIZE; WW++){

        double div_[3]={0,0,0};
        double w = 2.0*(WW-parm_.W_SIZE/2)*parm_.W_MAX/parm_.W_SIZE;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double kz = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
            for (int j = 0; j < parm_.K_SIZE; j++){
                double kx = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
                for (int l = 0; l < parm_.K_SIZE; l++){
                    double ky = 2.0*(l-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
                    if(fabs(kz)+fabs(kx)+fabs(ky)<=(3*pi/2)){
                        double k[3] = {kx,ky,kz};
                        Ham Ham_(parm_,k,0);
                        Opt_Green_transport_BI_div(parm_,Ham_,w,dk2,div_);
                    }
                    
                }
            }
        }

        for (int j = 0; j < 3; j++){
            div_sum[WW][j] = dw * FD(w,parm_.T)* div_[j] / (2*pi);
        }
    }

    
    
    for (int i = 0; i < SIZE; i++){
        for (int j = 0; j < 3; j++){
            div_total[j] += div_sum[i][j];
        }
        double w = 2.0*(i-parm_.W_SIZE/2)*parm_.W_MAX/parm_.W_SIZE;
        out1 << w;
        for (int j = 0; j < 3; j++){
            out1 << "\t" << div_sum[i][j];
        }
        out1 << endl;
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "div(x,y,z) " << div_total[0] << " " << div_total[1] << " " << div_total[2] << endl;
    //cout << "total " << Drude_sum+BCD_sum+Inj_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << parm_.delta;
    for (int i = 0; i < 3; i++){
        out3  << "\t" << div_total[i];
    }
    out3 << endl;
    return 0;
}