#include "transport.hpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <chrono>
#include <omp.h>

using namespace std;
typedef std::complex<double> Complex;

int main(int argc, char* argv[]){

    //ofstream out1("IFS_k.dat");
    ofstream out3("data.dat");

    parm parm_(argv);
    parm_.Parm_List();
    
    double dw = 2.0*parm_.W_MAX/parm_.W_SIZE;
    double dk2 = 4.0 * pi *pi / (parm_.K_SIZE*parm_.K_SIZE);
    
    cout << "dw " << dw << " dk " << dk2 <<endl;

    double gyro_xxy=0,gyro_yxy=0;
    double gyro_xxy_G=0,gyro_xxy_G1=0,gyro_xxy_G2=0,gyro_xxy_GWW=0,gyro_yxy_G=0,gyro_xxy_Gdiv=0,gyro_yxy_Gdiv=0;

    double gyro_GBI_xxy=0,gyro_GBIww_xxy=0;

    double gyro_xxy_sum[parm_.K_SIZE*parm_.K_SIZE],gyro_yxy_sum[parm_.K_SIZE*parm_.K_SIZE];
    //double IFS_xxy_G_sum[parm_.K_SIZE*parm_.K_SIZE],IFS_yxy_G_sum[parm_.K_SIZE*parm_.K_SIZE];
    //double gyro_xxy_G_sum[parm_.W_SIZE],gyro_yxy_G_sum[parm_.W_SIZE],gyro_xxy_Gdiv_sum[parm_.W_SIZE],gyro_yxy_Gdiv_sum[parm_.W_SIZE];
    double gyro_xxy_G_sum[parm_.W_SIZE],gyro_xxy_G1_sum[parm_.W_SIZE],gyro_xxy_G2_sum[parm_.W_SIZE],gyro_xxy_GWW_sum[parm_.W_SIZE];
    double gyro_xxy_GBI_sum[parm_.W_SIZE], gyro_xxy_GBIww_sum[parm_.W_SIZE];
    double im[M],re[M];
    im[0]=0;im[1]=0;re[0]=0;re[1]=0;

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    /*
    #pragma omp parallel for num_threads(52)
    for (int i = 0; i < parm_.K_SIZE; i++){
        double kx = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;

        for (int j = 0; j < parm_.K_SIZE; j++){
            double ky = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
            double k[D] = {kx,ky};
            Ham Ham_(parm_,k,0);
            //Opt_Gyration_BI(parm_,Ham_,gyro_xxy_sum[i+M*j],gyro_yxy_sum[i+M*j]);
            Opt_Gyration_BI3(parm_,Ham_,gyro_xxy_sum[i+parm_.K_SIZE*j],gyro_yxy_sum[i+parm_.K_SIZE*j]);
            //Opt_IFS_Green(parm_,Ham_,IFS_xxy_G_sum[i+M*j],IFS_yxy_G_sum[i+M*j]);
        }
    }

    
    
    for (int i = 0; i < parm_.K_SIZE; i++){
        for (int j = 0; j < parm_.K_SIZE; j++){
            //out1 << kx << "\t" << ky << "\t" << IFS_xxy_sum[i+M*j]/(4*pi*pi) << "\t" << IFS_yxy_sum[i+M*j]/(4*pi*pi)
            // << "\t" << IFS_xxy_G_sum[i+M*j]/(4*pi*pi) << "\t" << IFS_yxy_G_sum[i+M*j]/(4*pi*pi) << endl;
            gyro_xxy += dk2 * gyro_xxy_sum[i+parm_.K_SIZE*j] / (4*pi*pi);
            gyro_yxy += dk2 * gyro_yxy_sum[i+M*j] / (4*pi*pi);
            //IFS_xxy_G += dk2 * IFS_xxy_G_sum[i+M*j] / (4*pi*pi);
            //IFS_yxy_G += dk2 * IFS_yxy_G_sum[i+M*j] / (4*pi*pi);
        }
        //out1 << endl;
    }*/



    #pragma omp parallel for num_threads(52)
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        double w = (2*ww - parm_.W_SIZE)*parm_.W_MAX/parm_.W_SIZE;
        gyro_xxy_GBIww_sum[ww]=0;gyro_xxy_GBI_sum[ww]=0;
        gyro_xxy_G1_sum[ww]=0,gyro_xxy_G2_sum[ww]=0,gyro_xxy_GWW_sum[ww]=0,gyro_xxy_G_sum[ww]=0;
        for (int i = 0; i < parm_.K_SIZE; i++){
            double kx = 2.0*(i-parm_.K_SIZE/2)*pi/parm_.K_SIZE;

            for (int j = 0; j < parm_.K_SIZE; j++){
                double ky = 2.0*(j-parm_.K_SIZE/2)*pi/parm_.K_SIZE;
                double k[D] = {kx,ky};
                Ham Ham_(parm_,k,0);
                Opt_Circular_Green_BI(parm_,Ham_,w,gyro_xxy_GBIww_sum[ww],gyro_xxy_GBI_sum[ww]);
                Green Green_(parm_,w,im,re,Ham_.H_k);
                //Opt_Gyration_Green(parm_,Ham_,w,gyro_xxy_G_sum[ww],gyro_yxy_G_sum[ww]);
                Opt_Circular_Green(parm_,Ham_,Green_,w,gyro_xxy_G1_sum[ww],gyro_xxy_G2_sum[ww],gyro_xxy_GWW_sum[ww],gyro_xxy_G_sum[ww]);
            }
        }
    }
    for (int ww = 0; ww < parm_.W_SIZE; ww++){
        //gyro_xxy_G += dk2 * dw * gyro_xxy_G_sum[ww] / (4*pi*pi);
        //gyro_yxy_G += dk2 * dw * gyro_yxy_G_sum[ww] / (4*pi*pi);
        //IFS_xxy_Gdiv += dk2 * dw * IFS_xxy_Gdiv_sum[ww] / (4*pi*pi);
        //IFS_yxy_Gdiv += dk2 * dw * IFS_yxy_Gdiv_sum[ww] / (4*pi*pi);
        gyro_xxy_G1 += dw *gyro_xxy_G1_sum[ww] / (8*pi*pi*pi);
        gyro_xxy_G2 += dw *gyro_xxy_G2_sum[ww] / (8*pi*pi*pi);
        gyro_xxy_G += dw *gyro_xxy_G_sum[ww] / (8*pi*pi*pi);
        gyro_xxy_GWW += dw *gyro_xxy_GWW_sum[ww] / (8*pi*pi*pi);
        gyro_GBI_xxy += dw *gyro_xxy_GBI_sum[ww] / (8*pi*pi*pi);
        gyro_GBIww_xxy += dw *gyro_xxy_GBIww_sum[ww] / (8*pi*pi*pi);
    }


    end = chrono::system_clock::now();
    double time = static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0);


    cout << "=============================" << endl;
    cout <<"calculating time " << time << "min. (" << time * 60 << " sec.) " << endl;
    //cout << "IFS " << gyro_xxy << " " << gyro_yxy << endl;
    cout << "IFS " << gyro_GBIww_xxy << " " << gyro_GBI_xxy << endl;
    //cout << "IFS_Green " << gyro_xxy_G << " " << gyro_yxy_G << endl;
    cout << "IFS_Green " << gyro_xxy_G1 << " " << gyro_xxy_G2 << "  " << gyro_xxy_GWW << " " << gyro_xxy_G << endl;
    //cout << "total " << IFS_xxy_G + IFS_xxy_Gdiv << " " << IFS_yxy_G + IFS_yxy_Gdiv << endl;
    cout << "=============================" << endl;

    out3 << parm_.W << "\t" << parm_.mu << "\t" << gyro_xxy << "\t" << gyro_yxy << "\t" << gyro_xxy_G << endl;

    return 0;
}
