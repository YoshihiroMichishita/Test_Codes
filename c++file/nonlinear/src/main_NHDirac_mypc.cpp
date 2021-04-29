#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]){

    ofstream out("Dirac_NHDisp.dat");
    ofstream out2("Mresolved_Cond.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    //double im[parm_.W_SIZE][Mf],re[parm_.W_SIZE][Mf],
    double w[parm_.W_SIZE];
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = pi * 4.0 / (parm_.K_SIZE) * parm_.K_MAX / (parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;
    
    double im[Mf],re[Mf];
    for (int i = 0; i < Mf; i++){
        im[i]=0;re[i]=0;
    }
    

    //double XX_sum=0,XXX_sum=0,XXX2_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double XX_sum=0,YX_sum=0,XXXr_sum=0,XXXr2_sum=0,XXXi_sum=0,XXXi2_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double XX_sum2=0,XXX_sum2=0,XXX2_sum2=0,YX2_sum2=0,YXXr_sum2=0,YXXr2_sum2=0,YXXi_sum2=0,YXXi2_sum2=0;
    double n=0;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
 
    //cout << HSL_.SIZE << endl; 
    //const double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    
    for (int i = 0; i < parm_.K_SIZE; i++){
        double kx=2.0*(i-parm_.K_SIZE/2)*parm_.K_MAX/(parm_.K_SIZE);
        for (int j = 0; j < parm_.K_SIZE; j++){
            double ky=2.0*(j-parm_.K_SIZE/2)*parm_.K_MAX/(parm_.K_SIZE);
            double k[D]={kx,ky};
            Ham Ham_(parm_,k,2);
            double XX=0,YX=0,XXXr=0,XXXr2=0,XXXi=0,XXXi2=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0;
            double XX_=0,XXX_=0,XXX2_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0;
            #pragma omp parallel for num_threads(52) reduction(+: XX,XXXr,XXXr2,XXXi,XXXi2,YXXr,YXXr2,YXXi,YXXi2)
            for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                DC_NLH_NRC(dw,w[WW],parm_.T,Ham_,Green_,XX,YX,XXXr,XXXr2,XXXi,XXXi2,YXXr,YXXr2,YXXi,YXXi2);
            }
            /*#pragma omp parallel for num_threads(50) reduction(+: XX,XXX,XXX2,YXXr,YXXr2,YXXi,YXXi2,XX_,XXX_,XXX2_,YXXr_,YXXr2_,YXXi_,YXXi2_)
            for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
                //Linear_transport_NHwithNRC(parm_,dw,w[WW],parm_.T,Ham_,XX,XXX,XXX2,YXXr,YXXr2,YXXi,YXXi2,XX_,XXX_,XXX2_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }*/
            out2 << kx << "\t" << ky << "\t" << Ham_.NH_fac[0] << "\t" << Ham_.NH_fac[1] << "\t" << XX << "\t" << YX << "\t" << XXXi << "\t" << XXXi2 << "\t" << YXXi << "\t" << YXXi2 << endl;

            XX_sum += dk2 * XX / (4*pi*pi);
            YX_sum += dk2 * YX / (4*pi*pi);
            XXXi_sum += dk2 * XXXi / (4*pi*pi);
            XXXr_sum += dk2 * XXXr / (4*pi*pi);
            XXXi2_sum += dk2 * XXXi2 / (4*pi*pi);
            XXXr2_sum += dk2 * XXXr2 / (4*pi*pi);
            YXXi_sum += dk2 * YXXi / (4*pi*pi);
            YXXr_sum += dk2 * YXXr / (4*pi*pi);
            YXXi2_sum += dk2 * YXXi2 / (4*pi*pi);
            YXXr2_sum += dk2 * YXXr2 / (4*pi*pi);
            /*
            XX_sum += dk2 * XX / (4*pi*pi);
            XXX_sum += dk2 * XXX / (4*pi*pi);
            XXX2_sum += dk2 * XXX2 / (4*pi*pi);
            YXXi_sum += dk2 * YXXi / (4*pi*pi);
            YXXr_sum += dk2 * YXXr / (4*pi*pi);
            YXXi2_sum += dk2 * YXXi2 / (4*pi*pi);
            YXXr2_sum += dk2 * YXXr2 / (4*pi*pi);
            XX_sum2 += dk2 * XX_ / (4*pi*pi);
            XXX_sum2 += dk2 * XXX_ / (4*pi*pi);
            XXX2_sum2 += dk2 * XXX2_ / (4*pi*pi);
            YXXi_sum2 += dk2 * YXXi_ / (4*pi*pi);
            YXXr_sum2 += dk2 * YXXr_ / (4*pi*pi);
            YXXi2_sum2 += dk2 * YXXi2_ / (4*pi*pi);
            YXXr2_sum2 += dk2 * YXXr2_ / (4*pi*pi); */
            out << kx << "\t" << ky << "\t" << real(Ham_.E_NH[0]) << "\t" << real(Ham_.E_NH[1]) << endl;
        }
        out << endl;
        out2 << endl;
    }
    /*
    #pragma omp parallel for num_threads(52) reduction(+: XX_sum,XXX_sum,YXXr_sum,YXXr2_sum,YXXi_sum,YXXi2_sum,XX_sum2,XXX_sum2,YXXr_sum2,YXXr2_sum2,YXXi_sum2,YXXi2_sum2)
    for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
        double XX=0,XXX=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0;
        double XX_=0,XXX_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double kx=2.0*(i-parm_.K_SIZE/2)*pi/(parm_.K_SIZE);
            for (int j = 0; j < parm_.K_SIZE; j++){
                double ky=2.0*(j-parm_.K_SIZE/2)*pi/(parm_.K_SIZE);
                double k[D]={kx,ky};
                Ham Ham_(parm_,k,2);
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NHwithNRC(parm_,dw,w[WW],parm_.T,Ham_,XX,XXX,YXXr,YXXr2,YXXi,YXXi2,XX_,XXX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }
        XX_sum += dk2 * XX / (4*pi*pi);
        XXX_sum += dk2 * XXX / (4*pi*pi);
        YXXi_sum += dk2 * YXXi / (4*pi*pi);
        YXXr_sum += dk2 * YXXr / (4*pi*pi);
        YXXi2_sum += dk2 * YXXi2 / (4*pi*pi);
        YXXr2_sum += dk2 * YXXr2 / (4*pi*pi);
        
        
        XX_sum2 += dk2 * XX_ / (4*pi*pi);
        XXX_sum2 += dk2 * XXX_ / (4*pi*pi);
        YXXi_sum2 += dk2 * YXXi_ / (4*pi*pi);
        YXXr_sum2 += dk2 * YXXr_ / (4*pi*pi);
        YXXi2_sum2 += dk2 * YXXi2_ / (4*pi*pi);
        YXXr2_sum2 += dk2 * YXXr2_ / (4*pi*pi); 
        

        if(WW%(parm_.W_SIZE/10)==0){
            for (int p = 0; p < WW/(parm_.W_SIZE/10); p++){
                cout << "*";
            }
            cout << endl;
        }
    }*/

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);

    //double ZZ = 1.0-(re[i05][0]-re[im05][0])/(w[i05]-w[im05]);
    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    //cout << "n " << n << endl;
    cout << "beta(im,re)  (" << parm_.im_b << ", " << parm_.re_b << ")" << endl;
    //cout << "Z " << parm_.alpha << endl;
    //cout << "W " << parm_.W << endl;
    cout << "XX " << XX_sum << "\t" << XX_sum2 << endl;
    cout << "YX " << YX_sum << "\t" << endl;
    cout << "XXX(i,r) " << XXXi_sum << "\t" << XXXr_sum << endl;
    //cout << "XXX " << XXX_sum << "\t" << XXX_sum2 << endl;
    //cout << "XXX2 " << XXX2_sum << "\t" << XXX2_sum2 << endl;
    cout << "YXX " << YXXr_sum << "\t" << YXXi_sum << "\t" << YXXr_sum2 << "\t" << YXXi_sum2 << endl;
    cout << "YXX2 " << YXXr2_sum << "\t" << YXXi2_sum << "\t" << YXXr2_sum2 << "\t" << YXXi2_sum2 << endl;
    cout << "=============================" << endl;

    cout << parm_.im_b << "\t" << XX_sum << "\t" << 2.0*(2.0*XXXi_sum+XXXi2_sum) << "\t" << 2.0*(2*YXXi_sum+YXXi2_sum)
     << "\t" << XX_sum2 << "\t" << 2.0*(2.0*XXX_sum2+XXX2_sum2) << "\t" << 2.0*(2*YXXi_sum2+YXXi2_sum2) << endl;
    
    //out3 << parm_.im_b << "\t" << parm_.T << "\t" << XX_sum << "\t" << YX_sum << "\t" << 2.0*XXXi_sum << "\t" << 2.0*(2.0*XXXi_sum+XXXi2_sum) << "\t" << 2.0*YXXi_sum << "\t" << "\t" << 2.0*(2.0*YXXi_sum+YXXi2_sum) << "\t" << XX_sum2 << "\t" << 2.0*XXX_sum2 << "\t" << 2.0*(2.0*XXX_sum2+XXX2_sum2) << "\t" << 2.0*YXXi_sum2 << "\t" << 2.0*(2.0*YXXi_sum2+YXXi2_sum2) << endl;

    out3 << parm_.im_b << "\t" << parm_.T << "\t" << XX_sum << "\t" << YX_sum << "\t" << 2.0*XXXi_sum << "\t" << 2.0*(2.0*XXXi_sum+XXXi2_sum) << "\t" << 2.0*YXXi_sum << "\t" << "\t" << 2.0*(2.0*YXXi_sum+YXXi2_sum) << endl;

    return 0;
};