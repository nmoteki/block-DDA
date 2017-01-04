#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "fftw3.h"
#include "MBT_utils.hpp"
#include <omp.h>
using namespace Eigen;
using namespace std;

MatrixXcd bl_bicggr_jacobi_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
// Block BiCGGR [Tadano etal 2009 JSIAM letters]
Eigen::initParallel();
int num_element_occupy= address.size();
int num_element_cuboid= n.prod();
int L= B.cols();

VectorXcd jpre= 1.0/DIAG_A.array();

MatrixXcd B_pre= MatrixXcd::Zero(B.rows(),B.cols());
for(int l=0; l<L; ++l){
    B_pre.col(l)=B.col(l).array()*jpre.array();
}

  double B_prenorm= B_pre.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B_pre;
  MatrixXcd P= R;
  //MatrixXcd V= A*R;
  //// ---------------- FFT accerelation of V= A*R ------------------
  MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
  MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
  MatrixXcd V= MatrixXcd::Zero(B.rows(),L);
   // mvp_fft (diagonal)
   #pragma omp parallel for
   for(int l=0; l < L; ++l){
      Q1.col(l)=DIAG_A.array()*R.col(l).array();
       // mvp_fft (non-diagonal)
       for(int m=0; m < num_element_occupy; ++m){
           int mm= address[m];
           P_hat.col(l).segment(f*mm,f)= R.col(l).segment(f*m,f);
       }
       P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
       for(int m=0; m < num_element_occupy; ++m){
           int mm= address[m];
           Q1.col(l).segment(f*m,f)+= P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
       }
      V.col(l)= Q1.col(l).array()*jpre.array();
   }
   ///---------------------------------------------------------------


  MatrixXcd W= V;
  MatrixXcd R0til= R; //MatrixXcd::Random(n,L);
  MatrixXcd R0til_H= R0til.adjoint();

  MatrixXcd Y= MatrixXcd::Zero(B.rows(),L);

  for(int k= 0; k < itermax; ++k){
    MatrixXcd alpha= (R0til_H*V).partialPivLu().solve(R0til_H*R);
    complex<double> qsi= (W.adjoint()*R).trace()/(W.adjoint()*W).trace();
    MatrixXcd S= P-qsi*V;
    MatrixXcd U= S*alpha;
    //MatrixXcd Y= A*U;
    //// ---------------- FFT accerelation of Y= A*U ------------------
    P_hat.setZero();
     // mvp_fft (diagonal)
     #pragma omp parallel for
     for(int l=0; l < L; ++l){
        Q1.col(l)=DIAG_A.array()*U.col(l).array();
         // mvp_fft (non-diagonal)
         for(int m=0; m < num_element_occupy; ++m){
             int mm= address[m];
             P_hat.col(l).segment(f*mm,f)= U.col(l).segment(f*m,f);
         }
         P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
         for(int m=0; m < num_element_occupy; ++m){
             int mm= address[m];
             Q1.col(l).segment(f*m,f)+= P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
         }
         Y.col(l)= Q1.col(l).array()*jpre.array();
     }
     ///---------------------------------------------------------------

    X= X+qsi*R+U;
    MatrixXcd Rnew= R-qsi*W-Y;
    double err= Rnew.norm()/B_prenorm;
    cout << "bl_bicggr_jacobi_mvp_fft: " << "iter= " << k << " relative err= " << err << endl;
    if(err < tol) break;
    //W= A*Rnew;
    //// ---------------- FFT accerelation of W= A*Rnew ------------------

    P_hat.setZero();

    // mvp_fft (diagonal)
    #pragma omp parallel for
    for(int l=0; l < L; ++l){
       Q1.col(l)=DIAG_A.array()*Rnew.col(l).array();
        // mvp_fft (non-diagonal)
        for(int m=0; m < num_element_occupy; ++m){
            int mm= address[m];
            P_hat.col(l).segment(f*mm,f)= Rnew.col(l).segment(f*m,f);
        }
        P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
        for(int m=0; m < num_element_occupy; ++m){
            int mm= address[m];
            Q1.col(l).segment(f*m,f)+= P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
        }
        W.col(l)= Q1.col(l).array()*jpre.array();
    }
     ///---------------------------------------------------------------

    MatrixXcd gamma = (R0til_H*R).partialPivLu().solve(R0til_H*Rnew/qsi);
    R=Rnew;
    P= R+U*gamma;
    V= W+Y*gamma;

  }
  return X;
}
