#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
#include "fftw3.h"
#include "bl_MBT_utils.hpp"
using namespace Eigen;
using namespace std;
using namespace std::chrono;

MatrixXcd bl_bicgstab_mvp_fft(const VectorXi& n, const int& f, const vector<int>& address, const VectorXcd& Au_til, const VectorXcd& DIAG_A, const MatrixXcd& B, const double& tol, const int& itermax, fftw_plan& plan_fwd, fftw_plan& plan_inv)
{
    Eigen::initParallel();
// Block BiCGSTAB [Tadano etal 2009 JSIAM letters]
  int num_element_occupy= address.size();
  int num_element_cuboid= n.prod();
  int L= B.cols();


  double Bnorm= B.norm();
  MatrixXcd X= MatrixXcd::Zero(B.rows(),B.cols()); // Initial guess of X (zeros)
  MatrixXcd R= B;
  MatrixXcd P= R;
  MatrixXcd R0til= R; //MatrixXcd::Random(B.rows(),B.cols());
  MatrixXcd R0til_H= R0til.adjoint();


 for(int k= 0; k < itermax; ++k){

     //auto start_iter= steady_clock::now();
     //// ---------------- FFT accerelation of V= A*P ------------------
     MatrixXcd Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
      // mvp_fft (diagonal)
      for(int i=0; i < L; ++i){
         Q1.col(i)=DIAG_A.array()*P.col(i).array();
      }
      // mvp_fft (non-diagonal)
      MatrixXcd P_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
      for(int m=0; m < num_element_occupy; ++m){
          int mm= address[m];
          P_hat.block(f*mm,0,f,L)= P.block(f*m,0,f,L);
      }
      MatrixXcd Q_hat= bl_MBT_fft_mvp(n, f, Au_til, P_hat, plan_fwd, plan_inv);
      MatrixXcd Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
      for(int m=0; m < num_element_occupy; ++m){
          int mm= address[m];
          Q2.block(f*m,0,f,L)= Q_hat.block(f*mm,0,f,L); // non-diagonal contribution
      }

      MatrixXcd V= Q1+Q2;
      ///---------------------------------------------------------------


      FullPivLU<MatrixXcd> lu(R0til_H*V);
      //PartialPivLU<MatrixXcd> lu(R0til_H*V);
      MatrixXcd alfa= lu.solve(R0til_H*R);
      MatrixXcd T= R-V*alfa;


      // MatrixXcd Z= A*T;


      //// --------------- FFT accerelation of Z= A*T --------------------
      Q1= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
      // mvp_fft (diagonal)
      for(int i=0; i < L; ++i){
          Q1.col(i)=DIAG_A.array()*T.col(i).array();
      }
      // mvp_fft (non-diagonal)
      MatrixXcd T_hat= MatrixXcd::Zero(f*num_element_cuboid,L);
      for(int m=0; m < num_element_occupy; ++m){
          int mm= address[m];
          T_hat.block(f*mm,0,f,L)= T.block(f*m,0,f,L);
      }
      Q_hat= bl_MBT_fft_mvp(n, f, Au_til, T_hat, plan_fwd, plan_inv);

      Q2= MatrixXcd::Zero(B.rows(),L); // diagonal contribution
      for(int m=0; m < num_element_occupy; ++m){
          int mm= address[m];
          Q2.block(f*m,0,f,L)= Q_hat.block(f*mm,0,f,L); // non-diagonal contribution
      }
      MatrixXcd Z= Q1+Q2;
      ///---------------------------------------------------------------


      complex<double> qsi= (Z.adjoint()*T).trace()/(Z.adjoint()*Z).trace();
      X= X+P*alfa+qsi*T;
      R= T-qsi*Z;
      double err= R.norm()/Bnorm;
      cout << "bl_bicgstab: " << "iter= " << k << " relative err= " << err << endl;
      if(err < tol) break;
      MatrixXcd beta= lu.solve(-R0til_H*Z);
      P= R+(P-qsi*V)*beta;

     // auto end_iter= steady_clock::now();
     //cout << "1 iter time " << duration_cast<milliseconds>((end_iter - start_iter)).count() << " millisecond" << endl;
  }

  return X;
}
