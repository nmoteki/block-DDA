#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "fftw3.h"
#include "bl_MBT_utils.hpp"
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicggr_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
// Block BiCGGR [Tadano etal 2009 JSIAM letters]
int num_element_occupy= address.size();
int num_element_cuboid= n.prod();
int L= B.cols();

  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B;
  MatrixXcd P= R;
  //MatrixXcd V= A*R;
  //// ---------------- FFT accerelation of V= A*R ------------------
  MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
   // mvp_fft (diagonal)
   for(int i=0; i < L; ++i){
      Q1.col(i)=DIAG_A.array()*R.col(i).array();
   }
   // mvp_fft (non-diagonal)
   MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
   for(int m=0; m < num_element_occupy; ++m){
       int mm= address[m];
       P_hat.block(f*mm,0,f,L)= R.block(f*m,0,f,L);
   }
   MatrixXcd Q_hat= bl_MBT_fft_mvp(n, f, Au_til, P_hat, plan_fwd, plan_inv);
   MatrixXcd Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
   for(int m=0; m < num_element_occupy; ++m){
       int mm= address[m];
       Q2.block(f*m,0,f,L)= Q_hat.block(f*mm,0,f,L); // non-diagonal contribution
   }
   MatrixXcd V= Q1+Q2;
   ///---------------------------------------------------------------


  MatrixXcd W= V;
  MatrixXcd R0til= R; //MatrixXcd::Random(n,L);
  MatrixXcd R0til_H= R0til.adjoint();
  for(int k= 0; k < itermax; ++k){
    MatrixXcd alfa= (R0til_H*V).fullPivLu().solve(R0til_H*R);
    complex<double> qsi= (W.adjoint()*R).trace()/(W.adjoint()*W).trace();
    MatrixXcd S= P-qsi*V;
    MatrixXcd U= S*alfa;
    //MatrixXcd Y= A*U;
    //// ---------------- FFT accerelation of Y= A*U ------------------
     Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     // mvp_fft (diagonal)
     for(int i=0; i < L; ++i){
        Q1.col(i)=DIAG_A.array()*U.col(i).array();
     }
     // mvp_fft (non-diagonal)
     P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
     for(int m=0; m < num_element_occupy; ++m){
         int mm= address[m];
         P_hat.block(f*mm,0,f,L)= U.block(f*m,0,f,L);
     }
     Q_hat= bl_MBT_fft_mvp(n, f, Au_til, P_hat, plan_fwd, plan_inv);
     Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     for(int m=0; m < num_element_occupy; ++m){
         int mm= address[m];
         Q2.block(f*m,0,f,L)= Q_hat.block(f*mm,0,f,L); // non-diagonal contribution
     }
     MatrixXcd Y= Q1+Q2;
     ///---------------------------------------------------------------


    X= X+qsi*R+U;
    MatrixXcd Rnew= R-qsi*W-Y;
    double err= Rnew.norm()/Bnorm;
    cout << "bl_bicggr: " << "iter= " << k << " relative err= " << err << endl;
    if(err < tol) break;
    //W= A*Rnew;
    //// ---------------- FFT accerelation of W= A*Rnew ------------------
     Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     // mvp_fft (diagonal)
     for(int i=0; i < L; ++i){
        Q1.col(i)=DIAG_A.array()*Rnew.col(i).array();
     }
     // mvp_fft (non-diagonal)
     P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
     for(int m=0; m < num_element_occupy; ++m){
         int mm= address[m];
         P_hat.block(f*mm,0,f,L)= Rnew.block(f*m,0,f,L);
     }
     Q_hat= bl_MBT_fft_mvp(n, f, Au_til, P_hat, plan_fwd, plan_inv);
     Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     for(int m=0; m < num_element_occupy; ++m){
         int mm= address[m];
         Q2.block(f*m,0,f,L)= Q_hat.block(f*mm,0,f,L); // non-diagonal contribution
     }
     W= Q1+Q2;
     ///---------------------------------------------------------------


    MatrixXcd gamma = (R0til_H*R).fullPivLu().solve(R0til_H*Rnew/qsi);
    R=Rnew;
    P= R+U*gamma;
    V= W+Y*gamma;

  }


  return X;
}
