#include <complex>
#include <Eigen/Dense>
#include "fftw3.h"
#include "MBT_utils.hpp"
using namespace Eigen;
using namespace std;

VectorXcd BT_fft(const VectorXi& n, const int& f, MatrixXi n_ind, const int& level, const double& lf, const complex<double>& k) //テスト済み
{
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     n_ind                2 x M rectangular matrix indicating the position of each element in the MBT matrix
     level                current BT level (level=1,...,M).
     lf                   physical length per lattice spacing (length scale factor)
     k                   wavenumber in medium
    ***************************************************************************************/
    /*************************** Output ***************************************************
    MBT_elem              f x f general matrix of MBT level=M in Eigen::MatrixXcd
    ***************************************************************************************/

    VectorXcd Au;
    if(level == (n.size()+1)){
        // terminate recursion
        MatrixXcd Au1(f,f);
        Au.resize(f*f);
        Au1= application_function(n,f,n_ind,level,lf,k);
        Au1= Au1.colwise().reverse().eval(); // flipud in place
        Au1= Au1.transpose().eval(); // transpose in place !!! this is ncessary because Eigen is col-major order by defalt
        Map<VectorXcd> Au(Au1.data(),Au1.size());
        return Au;
    }else{
        //cout << "-------" << endl;
        int lenAu= f*f*(2*n.segment(level-1, n.size()-(level-1)).array()-1).prod();
        Au.resize(lenAu);
        int this_n= n(level-1);
        int b_edge= f*f*(2*n.segment(level, n.size()-level).array()-1).prod();
        //---lower triangular and diagonal blocks
        for(int i= this_n; i > 0; --i){
            n_ind(1,level-1)= 1;
            n_ind(0,level-1)= i;
            VectorXcd blk= BT_fft(n,f,n_ind,level+1,lf,k);
            Au.segment(b_edge*(this_n-i),b_edge)= blk;
        }
        //---upper triangular blocks
        for(int i= 2; i <= this_n; ++i){
            n_ind(1,level-1)=i;
            VectorXcd blk1= BT_fft(n,f,n_ind,level+1,lf,k);
            Au.segment(b_edge*(i+this_n-2),b_edge)= blk1;
        }
        return Au;
    }
}
