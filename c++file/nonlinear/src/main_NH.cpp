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
    ofstream out4("TMD_Dispw.dat");

    parm parm_(argv);
    parm_.Parm_List();

    /*
    ifstream Sigma; 
    Sigma.open(argv[15]);
    int NRG_SIZE=0;
    int im05,i05;
    double S0;
    double im0[parm_.W_SIZE][Mf],re0[parm_.W_SIZE][Mf],w0[parm_.W_SIZE];
    for( int i =0; i < parm_.W_SIZE; i++){
        Sigma >> w0[i] >> im0[i][0] >> re0[i][0];
        im0[i][1] = im0[i][0]; re0[i][1] = re0[i][0];
        if(w0[i]>-0.005 && w0[i-1]<-0.005){
            im05 = i;
        }
        else if(w0[i]>0 && w0[i-1]<0){
            S0 = (re0[i][0]+ re0[i-1][0])/2;
        }
        else if(w0[i]>0.005 && w0[i-1]<0.005){
            i05 = i;
        }
        if(!Sigma.good()){
            cout << i << endl;
            cout << im05 <<"\t"<< i05 <<endl;
            NRG_SIZE = i;
            break;
        }
    }*/

    //double im[parm_.W_SIZE][Mf],re[parm_.W_SIZE][Mf],
    double w[parm_.W_SIZE];
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }
    
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
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
    }

    for (int i = 0; i < HSL_.SIZE; i+=5){
        for (int WW = 0; WW < parm_.W_SIZE; WW+=10){
            out4 << HSL_.QL[i] << "\t" << w[WW] << "\t" << Spectral(parm_,w[WW],HSL_.KL[i],re[WW],im[WW]) << endl;
        }
        out4 << endl;
    }*/
    

    double XX_sum=0,YX_sum=0,YX2_sum=0,YXXr_sum=0,YXXr2_sum=0,YXXi_sum=0,YXXi2_sum=0;
    double XX_sum2=0,YX_sum2=0,YX2_sum2=0,YXXr_sum2=0,YXXr2_sum2=0,YXXi_sum2=0,YXXi2_sum2=0;
    double n=0;
    
    //chrono::system_clock::time_point start, end;
    //start = chrono::system_clock::now();
 
    //cout << HSL_.SIZE << endl; 
    const double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    #pragma omp parallel for reduction(+: XX_sum,YX_sum,YXXr_sum,YXXr2_sum,YXXi_sum,YXXi2_sum,XX_sum2,YX_sum2,YXXr_sum2,YXXr2_sum2,YXXi_sum2,YXXi2_sum2)
    for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
        double XX=0,YX=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0;
        double XX_=0,YX_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);
            for (int j = 0; j < parm_.K_SIZE; j++){
                double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
                Ham Ham_(parm_,k,2);
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
            
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
                Ham Ham_(parm_,k,2);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }

            double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
            for (int j = 1; j < parm_.K_SIZE+1; j++){
                double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
                double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
                Ham Ham_(parm_,k,2);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }
        XX_sum += dk2 * 4.0 * XX / (3*pi*pi)/(2 *pi);
        YX_sum += dk2 * 4.0 * YX / (3*pi*pi)/(2 *pi);
        YXXi_sum += dk2 * 4.0 * YXXi / (3*pi*pi)/(2 *pi);
        YXXr_sum += dk2 * 4.0 * YXXr / (3*pi*pi)/(2 *pi);
        YXXi2_sum += dk2 * 4.0 * YXXi2 / (3*pi*pi)/(2 *pi);
        YXXr2_sum += dk2 * 4.0 * YXXr2 / (3*pi*pi)/(2 *pi);
        
        
        XX_sum2 += dk2 * 4.0 * XX_ / (3*pi*pi)/(2 *pi);
        YX_sum2 += dk2 * 4.0 * YX_ / (3*pi*pi)/(2 *pi);
        YXXi_sum2 += dk2 * 4.0 * YXXi_ / (3*pi*pi)/(2 *pi);
        YXXr_sum2 += dk2 * 4.0 * YXXr_ / (3*pi*pi)/(2 *pi);
        YXXi2_sum2 += dk2 * 4.0 * YXXi2_ / (3*pi*pi)/(2 *pi);
        YXXr2_sum2 += dk2 * 4.0 * YXXr2_ / (3*pi*pi)/(2 *pi); 
        

        if(WW%(parm_.W_SIZE/10)==0){
            for (int p = 0; p < WW/(parm_.W_SIZE/10); p++){
                cout << "*";
            }
            cout << endl;
        }
    }

    /*
    #pragma omp parallel for reduction(+: XX_sum,YX_sum,YXXr_sum,YXXr2_sum,YXXi_sum,YXXi2_sum,XX_sum2,YX_sum2,YXXr_sum2,YXXr2_sum2,YXXi_sum2,YXXi2_sum2)
    for (int i = 0; i < parm_.K_SIZE; i++){
        double XX=0,YX=0,YXXr=0,YXXr2=0,YXXi=0,YXXi2=0;
        double XX_=0,YX_=0,YXXr_=0,YXXr2_=0,YXXi_=0,YXXi2_=0;
        double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);
        for (int j = 0; j < parm_.K_SIZE; j++){
            double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
            Ham Ham_(parm_,k,2);
            for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }
        
        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
            Ham Ham_(parm_,k,2);
            for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }

        double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
            Ham Ham_(parm_,k,2);
            for (int WW = 0; WW < parm_.W_SIZE-1; WW++){
                //Green Green_(parm_,w[WW],dw,im,re,Ham_.H_k);
                Linear_transport_NH(dw,w[WW],parm_.T,Ham_,XX,YX,YXXr,YXXr2,YXXi,YXXi2,XX_,YX_,YXXr_,YXXr2_,YXXi_,YXXi2_);
            }
        }
        XX_sum += dk2 * 4.0 * XX / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YX_sum += dk2 * 4.0 * YX / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXi_sum += dk2 * 4.0 * YXXi / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXr_sum += dk2 * 4.0 * YXXr / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXi2_sum += dk2 * 4.0 * YXXi2 / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXr2_sum += dk2 * 4.0 * YXXr2 / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        
        
        XX_sum2 += dk2 * 4.0 * XX_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YX_sum2 += dk2 * 4.0 * YX_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXi_sum2 += dk2 * 4.0 * YXXi_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXr_sum2 += dk2 * 4.0 * YXXr_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXi2_sum2 += dk2 * 4.0 * YXXi2_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi);
        YXXr2_sum2 += dk2 * 4.0 * YXXr2_ / (3*parm_.K_SIZE*parm_.K_SIZE)/(2 *pi); 
        if(i%(parm_.K_SIZE/10)==0){
            for (int p = 0; p < i/(parm_.K_SIZE/10); p++){
                cout << "*";
            }
            cout << endl;
        }
    }*/

    //end = chrono::system_clock::now();
    //double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);
    
    
    //double ZZ = 1.0-(re[i05][0]-re[im05][0])/(w[i05]-w[im05]);
    cout << "=============================" << endl;
    //cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    
    cout << "n " << n << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    cout << "XX " << XX_sum << "\t" << XX_sum2 << endl;
    cout << "YX " << YX_sum << "\t" << YX_sum2 << endl;
    cout << "YXX " << YXXr_sum << "\t" << YXXi_sum << "\t" << YXXr_sum2 << "\t" << YXXi_sum2 << endl;
    cout << "YXX " << YXXr2_sum << "\t" << YXXi2_sum << "\t" << YXXr2_sum2 << "\t" << YXXi2_sum2 << endl;
    cout << "=============================" << endl;
    cout << parm_.mu << "\t" << XX_sum << "\t" << YX_sum << "\t" << 2*YXXi_sum+YXXi2_sum
     << "\t" << XX_sum2 << "\t" << YX_sum2 << "\t" << 2*YXXi_sum2+YXXi2_sum2 << endl;
    
    out3 << parm_.NH << "\t" << parm_.alpha << "\t" << XX_sum << "\t" << YXXi_sum << "\t"
     << "\t" << 2.0*YXXi_sum+YXXi2_sum << "\t" << XX_sum2 << "\t" << YXXi_sum2 << "\t" << 2.0*YXXi_sum2+YXXi2_sum2 << endl;

    return 0;
};

