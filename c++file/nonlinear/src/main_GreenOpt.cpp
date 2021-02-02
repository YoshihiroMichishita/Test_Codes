//#include "matrix_op.hpp"
//#include "Ham_TMD.hpp"
#include "C3v_HSL.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out("TMD_Disp.dat");
    ofstream out1("TMD_Cond.dat");
    ofstream out2("Mresolved_Cond.dat");
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
    

    //double XX_sum=0,YX_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double PVXr_sum=0,PVXi_sum=0,PVYXr_sum=0,PVYXi_sum=0;
    double div1_sum=0,div2_sum=0,div3_sum=0;
    double n=0;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();



    const double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    for (int WW = 0; WW < (parm_.W_MAX+10*parm_.T) *parm_.W_SIZE/(2*parm_.W_MAX); WW++){
        #pragma omp parallel for reduction(+: PVXr_sum,PVYXr_sum,PVYXi_sum,div1_sum,div2_sum,div3_sum)
        for (int i = 0; i < parm_.K_SIZE; i++){
            double PVXr=0,PVXi=0,PVYXr=0,PVYXi=0;
            double div1=0,div2=0,div3=0;
            double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);
            for (int j = 0; j < parm_.K_SIZE; j++){
                double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Opt_transport(parm_,dw,w[WW],Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
                Ham_.~Ham();
                Green_.~Green();
            }
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Opt_transport(parm_,dw,w[WW],Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
                Ham_.~Ham();
                Green_.~Green();
            }

            double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
                Ham Ham_(parm_,k,1);
                Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Opt_transport(parm_,dw,w[WW],Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
                Ham_.~Ham();
                Green_.~Green();
            }
            
            PVXr_sum += 4.0 *dk2 * PVXr / (3*pi*pi)/(2*pi); 
            PVYXr_sum += 4.0 *dk2 * PVYXr / (3*pi*pi)/(2*pi); 
            PVYXi_sum += 4.0 *dk2 * PVYXi / (3*pi*pi)/(2*pi);

            div1_sum += 4.0 *dk2 * div1 / (3*pi*pi)/(2*pi); 
            div2_sum += 4.0 *dk2 * div2 / (3*pi*pi)/(2*pi);
            div3_sum += 4.0 *dk2 * div3 / (3*pi*pi)/(2*pi);
            
        }
            
        if(WW%(parm_.W_SIZE/10)==0){
            for (int l = 0; l < WW/(parm_.W_SIZE/10); l++){
                cout << "#"; //Parcentage of the finished work.
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
    cout << "W " << parm_.W << endl;
    cout << "PVX,DIV " << PVYXr_sum << "\t" << PVXr_sum << "\t" << div2_sum << endl;
    cout << "PVYX, DIV/w^2, DIV/w " << PVYXi_sum << "\t" << div1_sum << "\t" << div3_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << PVYXr_sum << "\t" << PVYXi_sum << "\t" << div1_sum << "\t" << div2_sum << "\t" << div3_sum << endl;
    return 0;
};

