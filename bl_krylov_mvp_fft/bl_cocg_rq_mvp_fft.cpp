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

MatrixXcd bl_cocg_rq_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
// only for complex-symmetric matrix
// Gu et al 2016, arXiv, Block variants of COCG and COCR methods
// for solving complex symmetric linear systems with multiple right-hand sides
Eigen::initParallel();
int num_element_occupy= address.size();
int num_element_cuboid= n.prod();
int L= B.cols();

  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)

  MatrixXcd Q;
  MatrixXcd xi;
  tie(Q,xi)= qr_reduced(B);
  MatrixXcd S= Q;

  for(int k= 0; k < itermax; ++k){
      //MatrixXcd AS= A*S;
     MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     MatrixXcd Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
     MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
     MatrixXcd Q_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
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
              Q_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
              for(int m=0; m < num_element_occupy; ++m){
                  int mm= address[m];
                  Q2.col(l).segment(f*m,f)= Q_hat.col(l).segment(f*mm,f); // non-diagonal contribution
              }
          }
      MatrixXcd AS= Q1+Q2;
      ///---------------------------------------------------------------


      MatrixXcd alpha= (S.transpose()*AS).fullPivLu().solve(Q.transpose()*Q);
      X= X+S*alpha*xi;
      MatrixXcd Qnew;
      MatrixXcd tau;

      tie(Qnew,tau)= qr_reduced(Q-AS*alpha);
      xi= tau*xi;
      double err= xi.norm()/Bnorm;
      cout << "bl_cocg_rq: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= (Q.transpose()*Q).fullPivLu().solve(tau.transpose()*Qnew.transpose()*Qnew);
      Q= Qnew;
      S= Q+S*beta;
  }

  return X;
}
