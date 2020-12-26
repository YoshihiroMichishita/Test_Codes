//#include "matrix_op.hpp"
//#include "Ham_TMD.hpp"
#include "C3v_HSL.hpp"
#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
//#include <omp.h>

using namespace std;


int main(int argc, char* argv[]){

    ofstream out("TMD_Disp.dat");
    ofstream out1("TMD_Cond.dat");
    ofstream out2("Mresolved_Cond.dat");
    ofstream out3("data.dat");
    ofstream out4("TMD_DOS.dat");

    parm parm_(argv);
    parm_.Parm_List();

    ifstream Sigma; 
    Sigma.open(argv[15]);
    int NRG_SIZE=0;
    int im05,i05;
    double S0;
    double dw = 2*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = pi * 4.0 / (3*parm_.K_SIZE) * 2.0 * pi / (sqrt(3.0)*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;
    /*
    HSL HSL_(parm_);

    for (int i = 0; i < HSL_.SIZE; i++){
        out << HSL_.QL[i];
        for (int WW = 0; WW < M; WW++){
            out << "\t" << HSL_.E[i][WW];
        }
        out << endl;
    }*/
    

    double Drude_sum=0,BCD_sum=0,Inj_sum=0,DrudeL_sum=0;
    int count=0;
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;
    
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
 
    double a1[D]={1.0,0},a2[D]={-0.5,sqrt(3.0)/2},a3[D]={-0.5,-sqrt(3.0)/2};
    for (int i = 0; i < parm_.K_SIZE; i++){

        double DrudeL_=0,Drude_=0,BCD_=0,Inj_=0;
        double k1[D]; k1[0]=4.0*i*a1[0]*pi/(3*parm_.K_SIZE); k1[1]=4.0*i*a1[1]*pi/(3*parm_.K_SIZE);

        for (int j = 0; j < parm_.K_SIZE; j++){
            double k2[D]; k2[0]=4.0*j*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*j*a2[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k2[0]; k[1]=k1[1]+k2[1];
            Ham Ham_(parm_,k,0);
            Opt_RTA_transport_BI2(parm_,Ham_,DrudeL_,Drude_,BCD_,Inj_);
        }

        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k1[0]+k3[0]; k[1]=k1[1]+k3[1];
            Ham Ham_(parm_,k,0);
            Opt_RTA_transport_BI2(parm_,Ham_,DrudeL_,Drude_,BCD_,Inj_);
        }

        double k2[D]; k2[0]=4.0*(i+1)*a2[0]*pi/(3*parm_.K_SIZE); k2[1]=4.0*(i+1)*a2[1]*pi/(3*parm_.K_SIZE);
        for (int j = 1; j < parm_.K_SIZE+1; j++){
            double k3[D]; k3[0]=4.0*j*a3[0]*pi/(3*parm_.K_SIZE); k3[1]=4.0*j*a3[1]*pi/(3*parm_.K_SIZE);
            double k[D]; k[0]=k2[0]+k3[0]; k[1]=k2[1]+k3[1];
            Ham Ham_(parm_,k,0);
            Opt_RTA_transport_BI2(parm_,Ham_,DrudeL_,Drude_,BCD_,Inj_);
        }
        DrudeL_sum += 4.0* dk2 * DrudeL_ / (3*pi*pi);
        Drude_sum += 4.0* dk2 * Drude_ / (3*pi*pi);
        BCD_sum += 4.0* dk2 * BCD_ / (3*pi*pi);
        Inj_sum += 4.0* dk2 * Inj_ / (3*pi*pi);
        
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "n " << count << endl;
    cout << "Z " << parm_.alpha << endl;
    cout << "W " << parm_.W << endl;
    cout << "DrudeL " << DrudeL_sum << endl;
    cout << "Drude " << Drude_sum << endl;
    cout << "BCD " << BCD_sum << endl;
    cout << "Inj " << Inj_sum << endl;
    cout << "total " << Drude_sum+BCD_sum+Inj_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.alpha << "\t" << parm_.W << "\t" << DrudeL_sum << "\t" << Drude_sum << "\t" << BCD_sum << "\t" << Inj_sum << endl;
    return 0;
};

