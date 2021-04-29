#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;
typedef std::complex<double> Complex;

int main(int argc, char* argv[]){

    ofstream out1("IFS_k.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    double dw = 2.0*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = 4.0 * pi *pi / (parm_.K_SIZE*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;

    double IFS_xxy=0,IFS_yxy=0;

    double IFS_xxy_sum[parm_.K_SIZE*parm_.K_SIZE],IFS_yxy_sum[parm_.K_SIZE*parm_.K_SIZE];
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    #pragma omp parallel for num_threads(52)
    for (int i = 0; i < parm_.K_SIZE; i++){

        double IFS_xxy=0,IFS_yxy=0;
        double kx = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;

        for (int j = 0; j < parm_.K_SIZE; j++){
            double ky = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
            double k[D] = {kx,ky};
            Ham Ham_(parm_,k,0);
            Opt_IFS_BI(parm_,Ham_,IFS_xxy_sum[i+M*j],IFS_yxy_sum[i+M*j]);
        }
    }

    
    
    for (int i = 0; i < parm_.K_SIZE; i++){
        double kx = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
        for (int j = 0; j < parm_.K_SIZE; j++){
            double ky = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
            out1 << kx << "\t" << ky << "\t" << IFS_xxy_sum[i+M*j]/(4*pi*pi) << "\t" << IFS_yxy_sum[i+M*j]/(4*pi*pi) << endl;
            IFS_xxy += dk2 * IFS_xxy_sum[i+M*j] / (4*pi*pi);
            IFS_yxy += dk2 * IFS_yxy_sum[i+M*j] / (4*pi*pi);
        }
        out1 << endl;
    }
    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    cout << "IFS " << IFS_xxy << " " << IFS_yxy << endl;
    //cout << "total " << Drude_sum+BCD_sum+Inj_sum << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << parm_.mu << "\t" << IFS_xxy << "\t" << IFS_yxy <<endl;;

    return 0;
}