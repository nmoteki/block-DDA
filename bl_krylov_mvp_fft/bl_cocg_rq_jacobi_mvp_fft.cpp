#include <iostream>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
#include "fftw3.h"
#include "MBT_utils.hpp"
#include <omp.h>
using namespace Eigen;
using namespace std;

MatrixXcd bl_cocg_rq_jacobi_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
// only for complex-symmetric matrix
// Gu et al 2016, arXiv, Block variants of COCG and COCR methods
// for solving complex symmetric linear systems with multiple right-hand sides
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

  MatrixXcd Q;
  MatrixXcd xi;
  tie(Q,xi)= qr_reduced(B_pre);
  MatrixXcd S= Q;

  MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
  MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
  MatrixXcd AS= MatrixXcd::Zero(B.rows(),L);

  for(int k= 0; k < itermax; ++k){
    
     P_hat.setZero();
     //// ---------------- FFT accerelation of AS= A*S ------------------
     // mvp_fft (diagonal)
     #pragma omp parallel for
         for(int l=0; l < L; ++l){
             Q1.col(l)=DIAG_A.array()*S.col(l).array();
              // mvp_fft (non-diagonal)
              for(int m=0; m < num_element_occupy; ++m){
                  int mm= address[m];
                  P_hat.col(l).segment(f*mm,f)= S.col(l).segment(f*m,f);
              }
              P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
              for(int m=0; m < num_element_occupy; ++m){
                  int mm= address[m];
                  Q1.col(l).segment(f*m,f)+= P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
              }
              AS.col(l)=Q1.col(l).array()*jpre.array();
          }
      ///---------------------------------------------------------------

      MatrixXcd alpha= (S.transpose()*AS).partialPivLu().solve(Q.transpose()*Q);
      X= X+S*alpha*xi;
      MatrixXcd Qnew;
      MatrixXcd tau;

      tie(Qnew,tau)= qr_reduced(Q-AS*alpha);
      xi= tau*xi;
      double err= xi.norm()/B_prenorm;
      cout << "bl_cocg_rq_jacobi_mvp_fft: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= (Q.transpose()*Q).partialPivLu().solve(tau.transpose()*Qnew.transpose()*Qnew);
      Q= Qnew;
      S= Q+S*beta;
  }

  return X;
}
