#include <tuple>
#include <complex>
# include <iostream>
#include <Eigen/Dense>
#include "linear_algebra_addon.hpp"
using namespace Eigen;
using namespace std;

tuple<MatrixXcd, MatrixXcd> qr_reduced(const MatrixXcd& A)
{
    // reduced QR decomposition using Modified Gram-Schmidt with Reortogonalization
    // reference: W. Ganter 1980, Algorithms for QR-Decomposition, research report no.80-02
    // A is  m x n  complex matrix
    int n= A.cols();
    MatrixXcd R= MatrixXcd::Zero(n,n);
    MatrixXcd Q= A;
    for(int k= 0; k < n; ++k){
        complex<double> tt {(0.0,0.0)};
        for(int j= 0; j < 2; ++j){
            for(int i= 0; i < k; ++i){
                complex<double> s= Q.col(i).adjoint()*Q.col(k);
                if(tt == (0.0,0.0)) R(i,k)= s;
                Q.col(k)= Q.col(k)-s*Q.col(i);
            }
            tt= Q.col(k).norm();
        }
        R(k,k)= tt;
        Q.col(k)= Q.col(k)/R(k,k);
    }
    return forward_as_tuple(Q,R);
}
