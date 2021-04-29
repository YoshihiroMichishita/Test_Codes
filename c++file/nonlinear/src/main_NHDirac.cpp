#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]){

    ofstream out("TMD_Disp.dat");
    ofstream out1("Mresolved_NH.dat");
    ofstream out2("Mresolved_Cond.dat");
    ofstream out3("data.dat");
    ofstream out4("TMD_Dispw.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    //double im[parm_.W_SIZE][Mf],re[parm_.W_SIZE][Mf],
    double w[parm_.W_SIZE];
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = pi * 4.0 / (parm_.K_SIZE) * pi / (parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;
    
    double im[Mf],re[Mf];
    for (int i = 0; i < Mf; i++){
        im[i]=0;re[i]=0;
    }
    

    double XX_sum=0,XXX_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double XX_sum2=0,XXX_sum2=0,YX2_sum2=0,YXXr_sum2=0,YXXr2_sum2=0,YXXi_sum2=0,YXXi2_sum2=0;
    double n=0;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
 
    //cout << HSL_.SIZE << endl; 
    //const double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    #pragma omp parallel for reduction(+: XX_sum,XXX_sum,YXXr_sum,YXXr2_sum,YXXi_sum,YXXi2_sum,XX_sum2,XXX_sum2,YXXr_sum2,YXXr2_sum2,YXXi_sum2,YXXi2_sum2)
    for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
        double XX=0,XXX=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0;
        double XX_=0,XXX_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double kx=2.0*(i-parm_.K_SIZE/2)*pi/(parm_.K_SIZE);
            for (int j = 0; j < parm_.K_SIZE; j++){
                double ky=2.0*(i-parm_.K_SIZE/2)*pi/(parm_.K_SIZE);
                double k[D]={kx,ky};
                Ham Ham_(parm_,k,2);
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NHwithNRC(dw,w[WW],parm_.T,Ham_,XX,XXX,YXXr,YXXr2,YXXi,YXXi2,XX_,XXX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }
        XX_sum += dk2 * XX / (4*pi*pi);
        XXX_sum += dk2 * XXX / (4*pi*pi);
        YXXi_sum += dk2 * YXXi / (4*pi*pi);
        YXXr_sum += dk2 * YXXr / (4*pi*pi);
        YXXi2_sum += dk2 * YXXi2 / (4*pi*pi);
        YXXr2_sum += dk2 * YXXr2 / (4*pi*pi);
        
        
        XX_sum2 += dk2 * XX_ / (4*pi*pi)/(2 *pi);
        XXX_sum2 += dk2 * XXX_ / (4*pi*pi)/(2 *pi);
        YXXi_sum2 += dk2 * YXXi_ / (4*pi*pi)/(2 *pi);
        YXXr_sum2 += dk2 * YXXr_ / (4*pi*pi)/(2 *pi);
        YXXi2_sum2 += dk2 * YXXi2_ / (4*pi*pi)/(2 *pi);
        YXXr2_sum2 += dk2 * YXXr2_ / (4*pi*pi)/(2 *pi); 
        

        if(WW%(parm_.W_SIZE/10)==0){
            for (int p = 0; p < WW/(parm_.W_SIZE/10); p++){
                cout << "*";
            }
            cout << endl;
        }
    }

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);
    
    
    //double ZZ = 1.0-(re[i05][0]-re[im05][0])/(w[i05]-w[im05]);
    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "n " << n << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    cout << "XX " << XX_sum << "\t" << XX_sum2 << endl;
    cout << "YX " << XXX_sum << "\t" << XXX_sum2 << endl;
    cout << "YXX " << YXXr_sum << "\t" << YXXi_sum << "\t" << YXXr_sum2 << "\t" << YXXi_sum2 << endl;
    cout << "YXX " << YXXr2_sum << "\t" << YXXi2_sum << "\t" << YXXr2_sum2 << "\t" << YXXi2_sum2 << endl;
    cout << "=============================" << endl;

    cout << parm_.mu << "\t" << XX_sum << "\t" << XXX_sum << "\t" << 2*YXXi_sum+YXXi2_sum
     << "\t" << XX_sum2 << "\t" << XXX_sum2 << "\t" << 2*YXXi_sum2+YXXi2_sum2 << endl;
    
    out3 << parm_.NH << "\t" << parm_.alpha << "\t" << XX_sum << "\t" << XXX_sum << "\t" << YXXi_sum << "\t"
     << "\t" << 2.0*YXXi_sum+YXXi2_sum << "\t" << XX_sum2 << "\t" << XXX_sum2 << "\t" << YXXi_sum2 << "\t" << 2.0*YXXi_sum2+YXXi2_sum2 << endl;

    return 0;
};