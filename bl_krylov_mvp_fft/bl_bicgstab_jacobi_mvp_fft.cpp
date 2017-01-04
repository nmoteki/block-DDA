#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
#include "fftw3.h"
#include "MBT_utils.hpp"
#include <omp.h>
using namespace Eigen;
using namespace std;
using namespace std::chrono;

MatrixXcd bl_bicgstab_jacobi_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
    Eigen::initParallel();
// Block BiCGSTAB [Tadano etal 2009 JSIAM letters]
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
  MatrixXcd R0til= R; //MatrixXcd::Random(B.rows(),B.cols());
  MatrixXcd R0til_H= R0til.adjoint();

  MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L);
  MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
  MatrixXcd V= MatrixXcd::Zero(B.rows(),L);
  MatrixXcd Z= MatrixXcd::Zero(B.rows(),L);

 for(int k= 0; k < itermax; ++k){

     //MatrixXcd V= A*P
     //// ---------------- FFT accerelation of V= A*P ------------------
     P_hat.setZero();
      // mvp_fft (diagonal)
      #pragma omp parallel for
      for(int l=0; l < L; ++l){
         Q1.col(l)=DIAG_A.array()*P.col(l).array();
          // mvp_fft (non-diagonal)
          for(int m=0; m < num_element_occupy; ++m){
              int mm= address[m];
              P_hat.col(l).segment(f*mm,f)= P.col(l).segment(f*m,f);
          }
          P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv); // in place
          for(int m=0; m < num_element_occupy; ++m){
              int mm= address[m];
              Q1.col(l).segment(f*m,f) += P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
          }
        V.col(l)= Q1.col(l).array()*jpre.array();
      }
      ///---------------------------------------------------------------

      PartialPivLU<MatrixXcd> lu(R0til_H*V);
      MatrixXcd alpha= lu.solve(R0til_H*R);
      MatrixXcd T= R-V*alpha;

      // MatrixXcd Z= A*T;
      //// --------------- FFT accerelation of Z= A*T --------------------
      P_hat.setZero();
      // mvp_fft (diagonal)
      #pragma omp parallel for
      for(int l=0; l < L; ++l){
         Q1.col(l)=DIAG_A.array()*T.col(l).array();
          // mvp_fft (non-diagonal)
          for(int m=0; m < num_element_occupy; ++m){
              int mm= address[m];
              P_hat.col(l).segment(f*mm,f)= T.col(l).segment(f*m,f);
          }
          P_hat.col(l)= MBT_fft_mvp(n, f, Au_til, P_hat.col(l), plan_fwd, plan_inv);
          for(int m=0; m < num_element_occupy; ++m){
              int mm= address[m];
              Q1.col(l).segment(f*m,f)+= P_hat.col(l).segment(f*mm,f); // non-diagonal contribution
          }
          Z.col(l)= Q1.col(l).array()*jpre.array();
      }
      ///---------------------------------------------------------------

      complex<double> qsi= (Z.adjoint()*T).trace()/(Z.adjoint()*Z).trace();
      X= X+P*alpha+qsi*T;
      R= T-qsi*Z;
      double err= R.norm()/B_prenorm;
      cout << "bl_bicgstab_jacobi_mvp_fft: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= lu.solve(-R0til_H*Z);
      P= R+(P-qsi*V)*beta;
  }
  return X;
}
