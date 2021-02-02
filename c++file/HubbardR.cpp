#include<cstdlib>
#include<fstream>
#include<vector>
#include<iostream>
#include<complex>
#include<math.h>
#include<string>
#include<sstream>
//#include<omp.h>
//#define MKL_Complex16 std::complex<double>
//#include<mkl.h>
//#include<mpi++.h>
//#define N 2

typedef std::complex<double> Complex;

using namespace std;



/*
void inverse(Complex H[N][N], Complex U[N][N]){
    int i,j;
    lapack_int lda=N;
    Complex a[lda*N];
    for( i=0; i<lda; i++ ){
        for( j=0; j<N; j++ ){
            a[lda*i+j] = H[i][j];
        }
    }
    lapack_int info1,info2;
    int ipiv[N];
    
    info1=LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, a, lda, ipiv);
    
    info2=LAPACKE_zgetri(LAPACK_ROW_MAJOR, N, a, lda, ipiv);
    
    for( i=0; i<lda; i++ ){
        for( j=0; j<N; j++ ){
            U[i][j]=a[lda*i+j];
        }
    }
    if (info2>0) {
        cout<<"error!"<<endl;
    }
}
*/


double mu_c;
double cutoff=1e-4;
double D;
//double DK;
double dk=M_PI*0.010;
//double a_ff;

double sign(double x){
  double result=0.0;
  if(x>0.0)
    result=1.0;
  if(x<0.0)
    result=-1.0;
  return result;
}


void Gii(double w,double im,double re,Complex& G, Complex& G0){
  //DOS between -1 and 1
  complex<double> I(0.0,1.0);
    //int info1,info2;
    //Complex H[N][N], U[N][N];
    //double E[N];
  //complex<double> x=(w-mu-re-I*im);
    G=0;
    for (double kx=0; kx<2*M_PI; kx+=dk) {
        for (double ky=0; ky<2*M_PI; ky+=dk) {
            double Ek=-2*D*(cos(kx)+cos(ky));
            G+=1.0/(w-(Ek+re+I*im+mu_c)+Delta*I);
        }
    }
    G0=1.0/G+re+im*I;
}





int main(int argc, char* argv[]){
  
  ofstream out1("Green_u.dat");
  ofstream out2("hyb_u.dat");
  ifstream Sigma_uu;
    
  //ifstream Sigma_u; //Sigmaからデータを取り込む
  //ifstream Sigma_d;

  if(argc<4){
    cout<<" Sigma  mu half-bandwidth cutoff alpha mu_f V"<<endl;
    exit(1);
  } //5個以上数値が入力されていない場合、警告を表示
  cout<<"open as uu "<<argv[1]<<endl;
    Sigma_uu.open(argv[1]);
    
  //Sigma_u.open(argv[1]); //Sigmaクラスでstart.impをopen
  //Sigma_d.open(argv[10]);

  
  string mystring=argv[2];
  stringstream stream(mystring); //streamにargv[2]の値を入れる
  stream>>mu_c;// muにargv[2]を代入
  cout<<"chemical potential of c "<<mu_c<<endl;//muは化学ポテンシャル
  mystring=argv[3];
  stringstream stream3(mystring);
  stream3>>D;//Dにargv[3]を代入
  cout<<"half bandwidth "<<D<<endl;//Dはバンド幅
  mystring=argv[4];
  stringstream stream2(mystring);
  stream2>>cutoff;//cutoffにargv[4]を代入
  cout<<"cutoff "<<cutoff<<endl;

  complex<double> I(0.0,1.0);

    //V=sqrt(0.01*J);
    
    int SIZE_W = 10000;
    
    double w[SIZE_W];
    double im[SIZE_W];
    double re[SIZE_W];
    int i;
    
    for (i=0; i<SIZE_W; i++) {
        Sigma_uu>>w[i]>>im[i]>>re[i];
        
        if(!Sigma_uu.good()){
            cout<<i<<endl;
            break;
        }
    }
    
    int SIZE = i;
    cout<<"size of w "<<SIZE<<endl;
    
    Complex g[SIZE];
    Complex h[SIZE];
    
    //#pragma omp parallel for
    for (int k=0; k<SIZE; k++) {
        if (fabs(w[k])>cutoff) {
            Gii(w[k],im[k],re[k],g[k],h[k]);
        }
    }
    
    for (int m=0; m<SIZE; m++) {
        if (fabs(w[m])>cutoff) {
            out1<<w[m]<<"\t"<<imag(g[m])<<"\t"<<real(g[m])<<endl;
            //out2<<w[m]<<"\t"<<fabs(imag(h[m]))<<"\t"<<real(h[m])<<endl;
            out2<<w[m]<<"\t"<<fabs(imag(h[m]))<<endl;
        }
  }
  return 0;  

}
