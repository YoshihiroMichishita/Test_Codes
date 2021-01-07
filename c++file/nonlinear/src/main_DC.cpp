//#include "matrix_op.hpp"
//#include "Ham_TMD.hpp"
#include "C3v_HSL.hpp"
#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;

int main(int argc, char* argv[]){

    ofstream out("TMD_Disp.dat");
    ofstream out1("TMD_Cond.dat");
    //ofstream out2("Mresolved_Cond.dat");
    ofstream out3("data.dat");
    ofstream out4("TMD_Dispw.dat");

    parm parm_(argv);
    parm_.Parm_List();

    ifstream Sigma; 
    Sigma.open(argv[15]);
    int NRG_SIZE=0;
    int im05,i05;
    double S0;

    double w[parm_.W_SIZE];
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }
    double dw = 2 * parm_.W_MAX/parm_.W_SIZE;
    double dk2 = pi * 4.0 / (3*parm_.K_SIZE) * 2.0 * pi / (sqrt(3.0)*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;
    
    double im[Mf],re[Mf];
    for (int i = 0; i < Mf; i++){
        im[i]=0;re[i]=0;
    }

    /*
    HSL HSL_(parm_);

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
            out4 << HSL_.QL[i] << "\t" << w[WW] << "\t" << Spectral(parm_,w[WW],HSL_.KL[i],re,im) << endl;
        }
        out4 << endl;
    }*/
    

    double XX_sum=0,YX_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double div1_sum=0,div2_sum=0,div3_sum=0;
    double n=0;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();


 
    //cout << HSL_.SIZE << endl; 

    double XX[parm_.W_SIZE],YX[parm_.W_SIZE],YXXr[parm_.W_SIZE],YXXi[parm_.W_SIZE],YXXr2[parm_.W_SIZE],YXXi2[parm_.W_SIZE];
    double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};

    #pragma omp parallel for reduction(+: XX_sum,YX_sum,YXXr_sum,YXXi_sum,YXXr2_sum,YXXi2_sum)
    for (int WW = 0; WW < parm_.W_SIZE; WW++){
        XX[WW]=0;YX[WW]=0;YXXr[WW]=0;YXXr2[WW]=0;YXXi[WW]=0;YXXi2[WW]=0;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);
            for (int j = 0; j < parm_.K_SIZE; j++){
                double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                DC_transport(dw,w[WW],parm_.T,Ham_,Green_,XX[WW],YX[WW],YXXr[WW],YXXr2[WW],YXXi[WW],YXXi2[WW]);
            }
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                DC_transport(dw,w[WW],parm_.T,Ham_,Green_,XX[WW],YX[WW],YXXr[WW],YXXr2[WW],YXXi[WW],YXXi2[WW]);
            }

            double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                DC_transport(dw,w[WW],parm_.T,Ham_,Green_,XX[WW],YX[WW],YXXr[WW],YXXr2[WW],YXXi[WW],YXXi2[WW]);
            }
            
        }
        XX_sum += 4.0 *dk2 * XX[WW] / (3*pi*pi)/(2*pi); 
        YX_sum += 4.0 *dk2 * YX[WW] / (3*pi*pi)/(2*pi);
        YXXi_sum += 4.0 *dk2 * YXXi[WW] / (3*pi*pi)/(2*pi); 
        YXXr_sum += 4.0 *dk2 * YXXr[WW] / (3*pi*pi)/(2*pi);
        YXXi2_sum += 4.0 *dk2 * YXXi2[WW] / (3*pi*pi)/(2*pi); 
        YXXr2_sum += 4.0 *dk2 * YXXr2[WW] / (3*pi*pi)/(2*pi);
        if(WW%(parm_.W_SIZE/10)==0){
            for (int l = 0; l < WW/(parm_.W_SIZE/10); l++){
                cout << "#";
            }
            cout << endl; 
        }
        
    }

    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);
    
    
    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "n " << n << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "XX " << XX_sum << endl;
    cout << "YX " << YX_sum << endl;
    cout << "YX " << YX2_sum << endl;
    cout << "YXX " << YXXr_sum << "\t" << YXXi_sum << endl;
    cout << "YXX2 " << YXXr2_sum << "\t" << YXXi2_sum << endl;
    cout << "=============================" << endl;
    out3 << parm_.mu + S0 << "\t" << n << "\t" << parm_.alpha << "\t" << XX_sum << "\t" << YXXi_sum << "\t" << YXXi2_sum << "\t" << 2*YXXi_sum+YXXi2_sum << endl;

    return 0;
};

