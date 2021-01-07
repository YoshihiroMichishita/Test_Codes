//#include "matrix_op.hpp"
//#include "Ham_TMD.hpp"
#include "C3v_HSL.hpp"
#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
//#include <chrono>
#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out("TMD_Disp.dat");
    ofstream out1("Mresolved_NH.dat");
    ofstream out2("Mresolved_Cond.dat");
    ofstream out3("data.dat");
    //ofstream out4("TMD_Dispw.dat");

    parm parm_(argv);
    parm_.Parm_List();

    ifstream Sigma; 
    Sigma.open(argv[15]);

    double w[parm_.W_SIZE];
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    double dk = 2*pi/parm_.K_SIZE;
    
    cout << "dw " << dw << " dk " << dk <<endl;
    
    double im[Mf],re[Mf];
    for (int i = 0; i < Mf; i++){
        im[i]=0;re[i]=0;
    }

    /*
    HSL HSL_(parm_);
    cout << HSL_.SIZE << endl;

    for (int i = 0; i < HSL_.SIZE; i++){
        out << HSL_.QL[i];
        for (int WW = 0; WW < M; WW++){
            out << "\t" << HSL_.E[i][WW];
        }
        out << endl;
    }*/
    /*
    for (int i = 0; i < HSL_.SIZE; i+=5){
        for (int WW = 0; WW < parm_.W_SIZE; WW+=10){
            out4 << HSL_.QL[i] << "\t" << w[WW] << "\t" << Spectral(parm_,w[WW],HSL_.KL[i],re[WW],im[WW]) << endl;
        }
        out4 << endl;
    }*/
    
    
    double XX_sum=0,YX_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0,XXX_sum=0;
    double XX_sum2=0,YX_sum2=0,YX2_sum2=0,YXXr_sum2=0,YXXr2_sum2=0,YXXi_sum2=0,YXXi2_sum2=0,XXX_sum2=0;

    double n=0;
    
    //chrono::system_clock::time_point start, end;
    //start = chrono::system_clock::now();
 
    //cout << HSL_.SIZE << endl; 
    for (int i = -parm_.K_SIZE/2; i < parm_.K_SIZE/2; i++){
        double kx = 2.0 * i * pi / (parm_.K_SIZE);
        double k[D]; 
        k[0] = kx; 
        //k[1] = ky;
        Ham Ham_(parm_,k,2);

        double XX=0,YX=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0,XXX=0;
        double XX_=0,YX_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0,XXX_=0;

        #pragma omp parallel for reduction(+: XX,XXX,XX_,XXX_)
        for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
            //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
            Linear_transport_NH_NRC(dw,w[WW],parm_.T,Ham_,XX,XXX,XX_,XXX_);
            
        }

        XX_sum += dk * XX /(2*pi)/(2*pi);
        XXX_sum += dk *  XXX/(2*pi)/(2*pi);
        XX_sum2 += dk *  XX_/(2*pi)/(2*pi);
        XXX_sum2 += dk * XXX_/(2*pi)/(2*pi);
        out2 << kx << "\t" << XX << "\t" << XX_ << "\t" << XXX << "\t" << XXX_ << endl; 
        out1 << kx << "\t" << Ham_.NH_fac[0] << "\t" << Ham_.NH_fac[1] << endl; 
        out1 << endl;
        out2 << endl;
    }
    //end = chrono::system_clock::now();
    //double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);
    
    
    //double ZZ = 1.0-(re[i05][0]-re[im05][0])/(w[i05]-w[im05]);
    cout << "=============================" << endl;
    //cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "n " << n << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    cout << "XX " << XX_sum << "\t" << XX_sum2 << endl;
    cout << "XXX " << XXX_sum << "\t" << XXX_sum2 << endl;
    cout << "=============================" << endl;
    cout << parm_.NH << "\t" << XX_sum << "\t" << XXX_sum << "\t" << XX_sum2 << "\t" << XXX_sum2 << endl;
    out3 << parm_.NH << "\t" << XX_sum << "\t" << XXX_sum << "\t" << XX_sum2 << "\t" << XXX_sum2 << endl;

    return 0;
};
