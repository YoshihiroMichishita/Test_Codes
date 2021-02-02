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


void k_sum_forNLC(parm parm_,double dk2, double w, double dw, double re[2], double im[2], double& PVXr,double& PVYXr,double& PVYXi,double& div1,double& div2,double& div3){
    #pragma omp parallel for reduction(+: PVXr,PVYXr,PVYXi,div1,div2,div3)
    for (int i = 0; i < parm_.K_SIZE; i++){
        double PVXr=0,PVXi=0,PVYXr=0,PVYXi=0;
        double div1=0,div2=0,div3=0;
        double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);
        for (int j = 0; j < parm_.K_SIZE; j++){
            double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
            Ham Ham_(parm_,k,1);
            Green Green_(parm_,w,dw,im,re,Ham_.H_k);
            Opt_transport(parm_,dw,w,Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
            Ham_.~Ham();
            Green_.~Green();
        }
        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
            Ham Ham_(parm_,k,1);
            Green Green_(parm_,w,dw,im,re,Ham_.H_k);
            Opt_transport(parm_,dw,w,Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
            Ham_.~Ham();
            Green_.~Green();
        }

        double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
            Ham Ham_(parm_,k,1);
            Green Green_(parm_,w,dw,im,re,Ham_.H_k);
            Opt_transport(parm_,dw,w,Ham_,Green_,PVYXr,PVXr,PVYXi,div1,div2,div3);
            Ham_.~Ham();
            Green_.~Green();
        }
        
        PVXr += 4.0 *dk2 * PVXr / (3*pi*pi)/(2*pi); 
        PVYXr += 4.0 *dk2 * PVYXr / (3*pi*pi)/(2*pi); 
        PVYXi += 4.0 *dk2 * PVYXi / (3*pi*pi)/(2*pi);

        div1 += 4.0 *dk2 * div1 / (3*pi*pi)/(2*pi); 
        div2 += 4.0 *dk2 * div2 / (3*pi*pi)/(2*pi);
        div3 += 4.0 *dk2 * div3 / (3*pi*pi)/(2*pi);
        
    }
    
};


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

    /*
    double w[parm_.W_SIZE];
    
    for (int i = 0; i < parm_.W_SIZE; i++){
        w[i] = 2.0*(i - parm_.W_SIZE/2) * parm_.W_MAX /(parm_.W_SIZE);
    }*/

    double dw = 2 * parm_.W_MAX/parm_.W_SIZE;
    double dk2 = pi * 4.0 / (3*parm_.K_SIZE) * 2.0 * pi / (sqrt(3.0)*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;
    
    double im[Mf],re[Mf];
    for (int i = 0; i < Mf; i++){
        im[i]=0;re[i]=0;
    }
    

    double PVXr_sum=0,PVXi_sum=0,PVYXr_sum=0,PVYXi_sum=0;
    double div1_sum=0,div2_sum=0,div3_sum=0;
    double n=0;
    double PVXr_old=0,PVXi_old=0,PVYXr_old=0,PVYXi_old=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    vector<double> w;
    vector<double> PVXr;
    vector<double> PVXi;
    vector<double> PVYXr;
    vector<double> PVYXi;
    double w_old=-parm_.W_MAX;

    const double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    for (double w0 = -parm_.W_MAX; w0 < parm_.W_MAX; w0 += dw){
        double PVXr_new=0,PVXi_new=0,PVYXr_new=0,PVYXi_new=0;
        w.push_back(w0);
        k_sum_forNLC(parm_,dk2,w0,dw,re,im,PVXr_new,PVYXr_new,PVYXi_new,div1_sum,div2_sum,div3_sum);
        PVXr.push_back(PVXr_new);
        PVXi.push_back(PVXi_new);
        PVYXr.push_back(PVYXr_new);
        PVYXi.push_back(PVYXi_new);
        PVXr_old=PVXr_new;PVXi_old=PVXi_new;
        PVYXr_old=PVYXr_new;PVYXi_old=PVYXi_new;
        while( fabs( PVXr_new - PVXr_old ) < parm_.delta/dw ){
            double w1 = (w.back()+w_old)/2;
            PVXr_new=0;PVXi_new=0;
            PVYXr_new=0;PVYXi_new=0;
            w.push_back(w1);
            k_sum_forNLC(parm_,dk2,w1,dw,re,im,PVXr_new,PVYXr_new,PVYXi_new,div1_sum,div2_sum,div3_sum);
            PVXr.insert(PVXr.end()-1,PVXr_new);
            PVXi.insert(PVXi.end()-1,PVXi_new);
            PVYXr.insert(PVYXr.end()-1,PVYXr_new);
            PVYXi.insert(PVYXi.end()-1,PVYXi_new);

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

