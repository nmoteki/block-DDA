#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
using namespace Eigen;
using namespace std;

int main(){
    FFT<double> fft;

    // fft for Eigen::Vector
    cout << " FFT in Eigen::Vector " << endl;
    VectorXcd in_Vec(10);
    for(int i= 0; i < 10; ++i) in_Vec(i)= complex<double>(1.0+0.0i)*static_cast<double>(i);

    VectorXcd out_Vec(10);
    cout << "in_Vec" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Vec(i) << endl;
    fft.fwd( out_Vec,in_Vec);
    cout << "out_Vec" << endl;
    for(int i= 0; i < 10; ++i) cout <<  out_Vec(i) << endl;
    fft.inv(in_Vec, out_Vec);
    cout << "in_Vec reconstruct" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Vec(i) << endl;

    cout << endl;

/*
    // fft for Eigen::Matrix
    cout << " FFT in Eigen::Matrix " << endl;
    MatrixXcd in_Mat(10,2);
    for(int i= 0; i < 10; ++i) in_Mat(i,1)= complex<double>(1.0+0.0i)*static_cast<double>(i);

    MatrixXcd out_Mat(10,2);
    VectorXcd tem_Mat(10);
    cout << "in_Mat" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Mat(i,1) << endl;
    fft.fwd(tem_Mat,in_Mat.col(1));
    out_Mat.col(1)= tem_Mat;
    cout << "out_Mat" << endl;
    for(int i= 0; i < 10; ++i) cout <<  out_Mat(i,1) << endl;
    fft.inv(tem_Mat, out_Mat.col(1));
    in_Mat.col(1)= tem_Mat;
    cout << "in_Mat reconstruct" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Mat(i,1) << endl;

    cout << endl;
*/

/*
    // fft for Eigen::Matrix
    cout << " FFT in Eigen::Matrix " << endl;
    MatrixXcd in_Mat(10,2);
    for(int i= 0; i < 10; ++i) in_Mat(i,1)= complex<double>(1.0+0.0i)*static_cast<double>(i);

    MatrixXcd out_Mat(10,2);
    cout << "in_Mat" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Mat(i,1) << endl;
    out_Mat.col(1)= fft.fwd(in_Mat.col(1));
    cout << "out_Mat" << endl;
    for(int i= 0; i < 10; ++i) cout <<  out_Mat(i,1) << endl;
    in_Mat.col(1)= fft.inv(out_Mat.col(1));
    cout << "in_Mat reconstruct" << endl;
    for(int i= 0; i < 10; ++i) cout << in_Mat(i,1) << endl;

    cout << endl;

*/

    // fft for Eigen::Matrix
    cout << " FFT in Eigen::Matrix " << endl;
    MatrixXcd in_Mat(10,2);
    for(int j= 0; j < 10; ++j) {
        in_Mat(j,0)= complex<double>(1.0+0.0i)*static_cast<double>(j);
        in_Mat(j,1)= complex<double>(2.0+0.0i)*static_cast<double>(j);
    }

    MatrixXcd out_Mat(10,2);
    cout << "in_Mat" << endl;
    cout << in_Mat << endl;
    for(int j=0; j < in_Mat.cols() ; ++j){
        out_Mat.col(j)= fft.fwd(in_Mat.col(j));
    }

    cout << "out_Mat" << endl;
    cout << out_Mat << endl;
    for(int j=0; j < out_Mat.cols(); ++j){
        in_Mat.col(j)= fft.inv(out_Mat.col(j));
    }

    cout << "in_Mat reconstruct" << endl;
    cout << in_Mat << endl;

    cout << endl;


    // fft for std::vector
    cout << " FFT in std::vector " << endl;
    vector<complex<double>> in_vec(10);
    for(int i= 0; i < 10; ++i) in_vec[i]= complex<double>(1.0+0.0i)*static_cast<double>(i);
    vector<complex<double>> out_vec(10,0);
    cout << "in_vec" << endl;
    for(int i= 0; i < 10; ++i) cout << in_vec[i] << endl;
    fft.fwd(out_vec,in_vec);
    cout << "out_vec" << endl;
    for(int i= 0; i < 10; ++i) cout << out_vec[i] << endl;
    fft.inv(in_vec,out_vec);
    cout << "in_vec reconstruct" << endl;
    for(int i= 0; i < 10; ++i) cout << in_vec[i] << endl;
}
