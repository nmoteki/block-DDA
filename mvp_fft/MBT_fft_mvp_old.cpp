#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include "MBT_utils.hpp"
using namespace Eigen;
using namespace std;

VectorXcd MBT_fft_mvp(const VectorXi& n,const int& f, const VectorXcd& Au_til, const VectorXcd& p_hat) //テスト済み
{
    //fast matrix-vector product using FFT. matrix is MBT
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     Au_til               fourier-transformed MBT projection of MBT matrix A
     p_hat                input vector
     /*************************** Output ***************************************************
     q_hat                result vector
     ***************************************************************************************/
     int fft_length= Au_til.size();
     //MBT Projection of data vector
     VectorXcd p_hatz= BT_pad(n, f, p_hat);
     //Zero padding to adjust xz size to fft_length
     int len_bt_pad= p_hatz.size();
     p_hatz.conservativeResize(fft_length);
     p_hatz.tail(fft_length-len_bt_pad).setZero();

     FFT<double> fft;
     VectorXcd p_hatz_til= fft.fwd(p_hatz);
     VectorXcd q_hatz_til= Au_til.cwiseProduct(p_hatz_til);
     VectorXcd q_hatz= fft.inv(q_hatz_til);
     VectorXcd q_hat= BT_reconstruct(n, f, q_hatz);

     return q_hat;
}
