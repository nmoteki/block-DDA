#include <iostream>
#include <complex>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

MatrixXcd application_function(const VectorXi& n, const int& f, MatrixXi& n_ind, int level, const double& lf, const complex<double>& k) //テスト済み
{
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     n_ind                2 x M rectangular matrix indicating the position of each element in the MBT matrix
     level                current BT level (level=1,...,M).
     lf                   physical length per lattice spacing (length scale factor)
     k                    wavenumber in medium
    ***************************************************************************************/
    /*************************** Output ***************************************************
    MBT_elem              f x f general matrix of MBT level=M in Eigen::MatrixXcd
    ***************************************************************************************/

    Vector3d Xkm= ((n_ind.row(0)-n_ind.row(1)).transpose()).cast<double>()*lf;
    double Xkm_abs= Xkm.norm();
    Matrix3d XkmXkm= Xkm*Xkm.transpose();


    switch (f)
    {
    case 3:
        {
            Matrix3cd MBT_elem;
            if(Xkm_abs == 0){
                MBT_elem= MatrixXcd::Zero(3,3);
                return MBT_elem;
            }
            Matrix3cd Gkm;
            Gkm= k*k*(Matrix3cd::Identity()-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k*Xkm_abs)/Xkm_abs/Xkm_abs)*(Matrix3cd::Identity()-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k*Xkm_abs)/Xkm_abs;
            MBT_elem= -Gkm;
            return MBT_elem;
        break;
        }
    case 6:
        {
            MatrixXcd MBT_elem(6,6);
            if(Xkm_abs == 0){
                MBT_elem= MatrixXcd::Zero(6,6);
                return MBT_elem;
            }
            Vector3d ukm= Xkm/Xkm_abs;
            Matrix3d ukm_vp= Matrix3d::Zero();
            ukm_vp(0,1)= -ukm(2);
            ukm_vp(0,2)= ukm(1);
            ukm_vp(1,0)= ukm(2);
            ukm_vp(1,2)= -ukm(0);
            ukm_vp(2,0)= -ukm(1);
            ukm_vp(2,1)= ukm(0);
            Matrix3cd Gkm;
            Gkm= k*k*(Matrix3cd::Identity()-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k*Xkm_abs)/Xkm_abs/Xkm_abs)*(Matrix3cd::Identity()-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k*Xkm_abs)/Xkm_abs;
            Matrix3cd Fkm;
            Fkm= -k*k*(exp(complex<double>(0.0+1.0i)*k*Xkm_abs)/Xkm_abs)*(1.0-1.0/(complex<double>(0.0+1.0i)*k*Xkm_abs))*ukm_vp;
            MBT_elem.block<3,3>(0,0)= -Gkm;
            MBT_elem.block<3,3>(0,3)= -Fkm;
            MBT_elem.block<3,3>(3,0)= Fkm;
            MBT_elem.block<3,3>(3,3)= -Gkm;
            return MBT_elem;
        break;
        }
    default:
        {
            cerr << "unknown f value :" << endl;
            exit(EXIT_FAILURE);
        break;
        }
    }
}// end of application_function()
