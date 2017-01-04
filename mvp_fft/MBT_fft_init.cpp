#include <iostream>
#include <complex>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include "fftw3.h"
#include "MBT_utils.hpp"
using namespace Eigen;
using namespace std;

VectorXcd MBT_fft_init(const VectorXi& n, const int& f, const double& lf, const complex<double>& k,fftw_plan& plan_fwd)
{
    //Prepare the fourier-transformed MBT projection of MBT matrix
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     lf                   physical length per lattice spacing (length scale factor)
     k                    wavenumber in medium
     /*************************** Output ***************************************************
     Au_til                result vector
     ***************************************************************************************/

    MatrixXi n_ind_init(2,n.size());
    n_ind_init= MatrixXi::Constant(2,n.size(),1);
    int fft_length= f*f*(2*n.array()-1).prod();

    int level_init= 1;
    //MBT projection of MBT matrix A;
    VectorXcd Au= BT_fft(n,f,n_ind_init,level_init,lf,k);
    if(!(fft_length == Au.size())){
        cerr << "wrong Au size !" << endl;
        exit(EXIT_FAILURE);
    }

    Au=Au.reverse().eval();
    fftw_execute_dft(plan_fwd,reinterpret_cast<fftw_complex*>(&Au[0]), reinterpret_cast<fftw_complex*>(&Au[0]));
    return Au;
    /*
    FFT<double> fft;
    VectorXcd Au_til= fft.fwd(Au);
    return Au_til;
    */
}
