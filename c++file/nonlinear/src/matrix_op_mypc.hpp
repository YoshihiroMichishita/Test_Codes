#pragma once

#include "const.hpp"
#include <iostream>

typedef std::complex<double> Complex;

using namespace std;

extern "C" {
    void zheev_(const char& JOBZ,const char& UPLO,
                const int& N, Complex** A, const int& LDA,
                double* W, Complex* WORK, const int& LWORK, double* RWORK,
                int& INFO,int JOBZlen,int UPLOlen );

    void zgeev_(const char& JOBVR,const char& JOBVL,
                const int& N, Complex* A, const int& LDA,
                Complex* W, Complex* VL, int& LDVL,Complex* VR, int& LDVR, Complex* WORK,
                const int& LWORK, double* RWORK,int& INFO);
    
    void zgetrf_( const int& S, const int& N, Complex** A,const int& LDA,
                 int* IPIV,int& INFO );
    
    void zgetri_( const int& N, Complex** A,const int& LDA,
                 int* IPIV, Complex* WORK, const int& LWORK,int& INFO );
};

//Matrix operating function

template <const int M> void inverse_NH(Complex H[M*M],Complex G[M*M]);
template <const int M> void Diag_H(Complex H[M*M], double E[M]);
template <const int M> void Diag_H(Complex H[M*M], Complex EV[M*M], double E[M]);
template <const int M> void Diag_NH(Complex H[M*M], Complex VL[M*M], Complex VR[M*M], Complex E[M]);
template <const int M> Complex Trace_NH(Complex H[M*M]);
template <const int M> double Trace_NHi(Complex H[M*M]);
template <const int M> double Trace_H(Complex H[M*M]);
template <const int M> double in(double A[M], double B[M]);
template <const int M> Complex InPro(Complex vL[M], Complex H[M*M], Complex vR[M]);
template <const int M> void InPro_M(Complex vL[M*M], Complex H[M*M], Complex vR[M*M],Complex H_n[M*M]);
template <const int M> void Prod2(Complex A1[M*M], Complex A2[M*M],Complex B[M*M]);
template <const int M> void Prod3(Complex A1[M*M], Complex A2[M*M], Complex A3[M*M],Complex B[M*M]);

template <const int M>
void inverse_NH(Complex H[M*M], Complex G[M*M]){

    int i,j;
    Complex U[M][M];
    for( i=0; i<M; i++ ){
        for( j=0; j<M; j++ ){
            U[i][j] = H[j*M+i];
        }
    }
    int info;
    int ipiv[4*M];
    const int lwork = 4*M;
    Complex work[4*M];
    
    zgetrf_( M, M, (Complex**)U, M, ipiv, info);
    
    zgetri_( M, (Complex**)U, M, ipiv, work, lwork, info);
    
    for( i=0; i<M; i++ ){
        for( j=0; j<M; j++ ){
            G[i*M+j] = U[j][i];
        }
    }
    
    if (info>0) {
        std::cout<<"error!"<<endl;
    }
}

template <const int M>
void Diag_H(Complex H[M*M], double E[M]){
    int i,j;
    int lda=M;
    Complex a[M][M];
    for( i=0; i<lda; i++ ){
        E[i]=0;
        for( j=0; j<M; j++ ){
            a[j][i] = H[i*M+j];
        }
    }
    int info1;
    Complex work[4*M];
    const int lwork=4*M;
    double rwork[4*M];
    zheev_('N','U',M,(Complex**)a,lda,E,work,lwork,rwork,info1,1,1);
    if (info1>0) {
        cout<<"error!"<<endl;
    }
}

template <const int M>
void Diag_H(Complex H[M*M],Complex EV[M*M], double E[M]){

    int i,j;
    int lda=M;
    Complex a[M][M];
    for( i=0; i<lda; i++ ){
        E[i]=0;
        for( j=0; j<M; j++ ){
            a[j][i] = H[i*M+j];
        }
    }
    int info1;
    Complex work[4*M];
    const int lwork=4*M;
    double rwork[4*M];
    zheev_('V','U',M,(Complex**)a,lda,E,work,lwork,rwork,info1,1,1);
    if (info1>0) {
        cout<<"error!"<<endl;
    }
    for( i=0; i<lda; i++ ){
        for( j=0; j<M; j++ ){
            EV[lda*i+j] = a[j][i];
        }
    }

}

template <const int M> 
void Diag_NH(Complex H[M*M], Complex VL[M*M], Complex VR[M*M], Complex E[M]){
    int i,j;
    int lda=M;
    int ldvl=M;
    int ldvr=M;
    
    int info1;
    
    Complex work[4*M];
    const int lwork=4*M;
    double rwork[4*M];

    Complex a[M*M];
    for( i=0; i<M; i++ ){
        for( j=0; j<M; j++ ){
            a[j*M+i] = H[i*M+j];
        }
    };
    Complex vl[M*M],vr[M*M];
    
    zgeev_('V','V', M, a, lda, E,vl,ldvl,vr,ldvr,work,lwork,rwork,info1);
    
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            VL[i*M+j] = conj(vl[i*M+j]);
            VR[i*M+j] = vr[i+j*M];
        }
    }

    if (info1>0) {
        std::cout<<"error!"<<endl;
    }
    if (info1<0) {
        std::cout<<info1<<"-th parameter has the illegal value" <<endl;
    }
};

template <const int M>
Complex Trace_NH(Complex H[M*M]){
    Complex C=0;
    for (int i = 0; i < M; i++){
        C += H[i*(M+1)];
    }
    return C;
}

template <const int M>
double Trace_NHi(Complex H[M*M]){
    double C=0;
    for (int i = 0; i < M; i++){
        C += imag(H[i*(M+1)]);
    }
    return C;
}

template <const int M>
double Trace_H(Complex H[M*M]){
    double C=0;
    for (int i = 0; i < M; i++){
        C += real(H[i*(M+1)]);
    }
    return C;
}

template <const int M> 
double in(double A[M], double B[M]){
    double IN = 0;
    for (int i = 0; i < M; i++){
        IN += A[i]*B[i];
    }
    return IN;
}

template <const int M>
Complex InPro(Complex vL[M], Complex H[M*M], Complex vR[M]){
    Complex C=0;
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            C += vL[i]*H[i*M+j]*vR[j];
        }
    }
    return C;
}

template <const int M>
void InPro_M(Complex vL[M*M], Complex H[M*M], Complex vR[M*M],Complex H_n[M*M]){
    for (int k = 0; k < M; k++){
        for (int l = 0; l < M; l++){
            H_n[k*M+l]=0;
            for (int i = 0; i < M; i++){
                for (int j = 0; j < M; j++){
                    H_n[k*M+l] += vL[k+M*i]*H[i*M+j]*vR[M*j+l];
                }
            }
        }     
    }
};

template <const int M>
void Prod2(Complex A1[M*M],Complex A2[M*M],Complex B[M*M]){
    for (int i = 0; i < M; i++){
        for (int j = 0; j < M; j++){
            B[i*M+j] = 0;
            for (int k = 0; k < M; k++){
                B[i*M+j] += A1[i*M+k]*A2[k*M+j];
            }
        }
    }
};

template <const int M>
void Prod3(Complex A1[M*M],Complex A2[M*M],Complex A3[M*M], Complex B[M*M]){

    for (int i = 0; i < M*M; i++){
        B[i] = 0;
    }
    Complex C[M*M];
    for (int i = 0; i < M*M; i++){
        C[i] = 0;
    }

    Prod2<M>(A1,A2,C);
    Prod2<M>(C,A3,B);

};
