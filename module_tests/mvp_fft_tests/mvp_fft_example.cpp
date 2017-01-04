#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <chrono>
//#include <ctime>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
using namespace Eigen;
using namespace std;
using namespace std::chrono;




MatrixXcd application_function(const VectorXi& n, const int& f, MatrixXi& n_ind, int level, const double& lf, const double& k0) //テスト済み
{
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     n_ind                2 x M rectangular matrix indicating the position of each element in the MBT matrix
     level                current BT level (level=1,...,M).
     lf                   physical length per lattice spacing (length scale factor)
     k0                   free space wavenumber
    ***************************************************************************************/
    /*************************** Output ***************************************************
    MBT_elem              f x f general matrix of MBT level=M in Eigen::MatrixXcd
    ***************************************************************************************/

    Vector3d Xkm= ((n_ind.row(0)-n_ind.row(1)).transpose()).cast<double>()*lf;
    double Xkm_abs= Xkm.norm();
    Matrix3d XkmXkm= Xkm*Xkm.transpose();


    switch (f)
    {
    case 1:
        {
            MatrixXcd MBT_elem(1,1);

            MatrixXi n_ind_222_111(2,3);  n_ind_222_111 << 2,2,2,1,1,1;
            MatrixXi n_ind_221_111(2,3);  n_ind_221_111 << 2,2,1,1,1,1;
            MatrixXi n_ind_221_112(2,3);  n_ind_221_112 << 2,2,1,1,1,2;
            MatrixXi n_ind_212_111(2,3);  n_ind_212_111 << 2,1,2,1,1,1;
            MatrixXi n_ind_211_111(2,3);  n_ind_211_111 << 2,1,1,1,1,1;
            MatrixXi n_ind_211_112(2,3);  n_ind_211_112 << 2,1,1,1,1,2;
            MatrixXi n_ind_212_121(2,3);  n_ind_212_121 << 2,1,2,1,2,1;
            MatrixXi n_ind_211_121(2,3);  n_ind_211_121 << 2,1,1,1,2,1;
            MatrixXi n_ind_211_122(2,3);  n_ind_211_122 << 2,1,1,1,2,2;
            MatrixXi n_ind_122_111(2,3);  n_ind_122_111 << 1,2,2,1,1,1;
            MatrixXi n_ind_121_111(2,3);  n_ind_121_111 << 1,2,1,1,1,1;
            MatrixXi n_ind_121_112(2,3);  n_ind_121_112 << 1,2,1,1,1,2;
            MatrixXi n_ind_112_111(2,3);  n_ind_112_111 << 1,1,2,1,1,1;
            MatrixXi n_ind_111_111(2,3);  n_ind_111_111 << 1,1,1,1,1,1;
            MatrixXi n_ind_111_112(2,3);  n_ind_111_112 << 1,1,1,1,1,2;
            MatrixXi n_ind_112_121(2,3);  n_ind_112_121 << 1,1,2,1,2,1;
            MatrixXi n_ind_111_121(2,3);  n_ind_111_121 << 1,1,1,1,2,1;
            MatrixXi n_ind_111_122(2,3);  n_ind_111_122 << 1,1,1,1,2,2;
            MatrixXi n_ind_122_211(2,3);  n_ind_122_211 << 1,2,2,2,1,1;
            MatrixXi n_ind_121_211(2,3);  n_ind_121_211 << 1,2,1,2,1,1;
            MatrixXi n_ind_121_212(2,3);  n_ind_121_212 << 1,2,1,2,1,2;
            MatrixXi n_ind_112_211(2,3);  n_ind_112_211 << 1,1,2,2,1,1;
            MatrixXi n_ind_111_211(2,3);  n_ind_111_211 << 1,1,1,2,1,1;
            MatrixXi n_ind_111_212(2,3);  n_ind_111_212 << 1,1,1,2,1,2;
            MatrixXi n_ind_112_221(2,3);  n_ind_112_221 << 1,1,2,2,2,1;
            MatrixXi n_ind_111_221(2,3);  n_ind_111_221 << 1,1,1,2,2,1;
            MatrixXi n_ind_111_222(2,3);  n_ind_111_222 << 1,1,1,2,2,2;

            if((n_ind-n_ind_222_111).norm()==0) MBT_elem <<1;
            if((n_ind-n_ind_221_111).norm()==0) MBT_elem <<2;
            if((n_ind-n_ind_221_112).norm()==0) MBT_elem <<3;
            if((n_ind-n_ind_212_111).norm()==0) MBT_elem <<4;
            if((n_ind-n_ind_211_111).norm()==0) MBT_elem <<5;
            if((n_ind-n_ind_211_112).norm()==0) MBT_elem <<6;
            if((n_ind-n_ind_212_121).norm()==0) MBT_elem <<7;
            if((n_ind-n_ind_211_121).norm()==0) MBT_elem <<8;
            if((n_ind-n_ind_211_122).norm()==0) MBT_elem <<9;
            if((n_ind-n_ind_122_111).norm()==0) MBT_elem <<10;
            if((n_ind-n_ind_121_111).norm()==0) MBT_elem <<11;
            if((n_ind-n_ind_121_112).norm()==0) MBT_elem <<12;
            if((n_ind-n_ind_112_111).norm()==0) MBT_elem <<13;
            if((n_ind-n_ind_111_111).norm()==0) MBT_elem <<14;
            if((n_ind-n_ind_111_112).norm()==0) MBT_elem <<15;
            if((n_ind-n_ind_112_121).norm()==0) MBT_elem <<16;
            if((n_ind-n_ind_111_121).norm()==0) MBT_elem <<17;
            if((n_ind-n_ind_111_122).norm()==0) MBT_elem <<18;
            if((n_ind-n_ind_122_211).norm()==0) MBT_elem <<19;
            if((n_ind-n_ind_121_211).norm()==0) MBT_elem <<20;
            if((n_ind-n_ind_121_212).norm()==0) MBT_elem <<21;
            if((n_ind-n_ind_112_211).norm()==0) MBT_elem <<22;
            if((n_ind-n_ind_111_211).norm()==0) MBT_elem <<23;
            if((n_ind-n_ind_111_212).norm()==0) MBT_elem <<24;
            if((n_ind-n_ind_112_221).norm()==0) MBT_elem <<25;
            if((n_ind-n_ind_111_221).norm()==0) MBT_elem <<26;
            if((n_ind-n_ind_111_222).norm()==0) MBT_elem <<27;

            //cout << n_ind << endl;
            //cout << "MBT_elem= " << MBT_elem << endl;

            return MBT_elem;

        break;
        }
    case 3:
        {
            Matrix3cd MBT_elem;
            if(Xkm_abs == 0){
                MBT_elem= MatrixXcd::Zero(3,3);
                return MBT_elem;
            }
            Matrix3cd Gkm;
            Gkm= k0*k0*(Matrix3cd::Identity()-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs/Xkm_abs)*(Matrix3cd::Identity()-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs;
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
            Gkm= k0*k0*(Matrix3cd::Identity()-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs/Xkm_abs)*(Matrix3cd::Identity()-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
            Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs;
            Matrix3cd Fkm;
            Fkm= -k0*k0*(exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs)*(1.0-1.0/(complex<double>(0.0+1.0i)*k0*Xkm_abs))*ukm_vp;
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



VectorXcd BT_fft(const VectorXi& n, const int& f, MatrixXi n_ind, const int& level, const double& lf, const double& k0) //テスト済み
{
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     n_ind                2 x M rectangular matrix indicating the position of each element in the MBT matrix
     level                current BT level (level=1,...,M).
     lf                   physical length per lattice spacing (length scale factor)
     k0                   free space wavenumber
    ***************************************************************************************/
    /*************************** Output ***************************************************
    MBT_elem              f x f general matrix of MBT level=M in Eigen::MatrixXcd
    ***************************************************************************************/

    VectorXcd Au;
    if(level == (n.size()+1)){
        // terminate recursion
        MatrixXcd Au1(f,f);
        Au.resize(f*f);
        Au1= application_function(n,f,n_ind,level,lf,k0);
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
            VectorXcd blk= BT_fft(n,f,n_ind,level+1,lf,k0);
            Au.segment(b_edge*(this_n-i),b_edge)= blk;
        }
        //---upper triangular blocks
        for(int i= 2; i <= this_n; ++i){
            n_ind(1,level-1)=i;
            VectorXcd blk1= BT_fft(n,f,n_ind,level+1,lf,k0);
            Au.segment(b_edge*(i+this_n-2),b_edge)= blk1;
        }
        return Au;
    }
}


VectorXcd BT_pad(const VectorXi& n, const int& f, const VectorXcd& x) // テスト済み
{
//Generation of xz by inserting zeros into x
/*************************** Input parameters *****************************************
 n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
 f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
 x                    input vector before padding
 /*************************** Output ***************************************************
 xz                   padded vector
 ***************************************************************************************/
    VectorXcd xz;
    int lenxz;
    if(n.size() < 1){
        lenxz= x.size()+f*f-f;
        xz.resize(lenxz);
        xz.segment(0,x.size())= x;
        xz.segment(x.size(),lenxz-x.size())= VectorXcd::Zero(lenxz-x.size());
        return xz;
    }else{
        int b_edge= f*(n.segment(1, n.size()-1)).prod();
        int this_n= n(0);

        for(int i= 1; i <= this_n; ++i){
            int lenxz_old= xz.size();
            int lenblk= BT_pad(n.segment(1, n.size()-1), f, x.segment((i-1)*b_edge,b_edge)).size();
            if(n.size() > 1){
                lenxz= lenxz_old+lenblk+(n(1)-1)*f*f*(2*n.segment(2,(n.size()-2)).array()-1).prod();
                xz.conservativeResize(lenxz);
                xz.segment(lenxz_old, xz.size()-lenxz_old)= VectorXcd::Zero(xz.size()-lenxz_old);
                xz.segment(lenxz_old,lenblk)= BT_pad(n.segment(1, n.size()-1), f, x.segment((i-1)*b_edge,b_edge));
            }else{
                lenxz= lenxz_old+lenblk;
                xz.conservativeResize(lenxz);
                xz.segment(lenxz_old, xz.size()-lenxz_old)= BT_pad(n.segment(1, n.size()-1), f, x.segment((i-1)*b_edge,b_edge));
            }
        }
        return xz;
    }

}


MatrixXcd bl_BT_pad(const VectorXi& n, const int& f, const MatrixXcd& x) // テスト済み
{
//Generation of xz by inserting zeros into x
/*************************** Input parameters *****************************************
 n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
 f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
 x                    input (n x L) matrix before padding
 /*************************** Output ***************************************************
 xz                   padded matrix
 ***************************************************************************************/
    MatrixXcd xz;
    int L= x.cols();
    int lenxz;
    if(n.size() < 1){
        lenxz= x.rows()+f*f-f;
        xz.resize(lenxz,L);
        xz.block(0,0,x.rows(),L)= x;
        xz.block(x.rows(),0,lenxz-x.rows(),L).setZero();
        return xz;
    }else{
        int b_edge= f*(n.segment(1, n.size()-1)).prod();
        int this_n= n(0);

        for(int i= 1; i <= this_n; ++i){
            int lenxz_old= xz.rows();
            int lenblk= bl_BT_pad(n.segment(1, n.size()-1), f, x.block((i-1)*b_edge, 0, b_edge, L)).rows();
            if(n.size() > 1){
                lenxz= lenxz_old+lenblk+(n(1)-1)*f*f*(2*n.segment(2,(n.size()-2)).array()-1).prod();
                xz.conservativeResize(lenxz,L);
                xz.block(lenxz_old, 0, xz.rows()-lenxz_old, L).setZero();
                xz.block(lenxz_old, 0, lenblk, L)= bl_BT_pad(n.segment(1, n.size()-1), f, x.block((i-1)*b_edge, 0, b_edge, L));
            }else{
                lenxz= lenxz_old+lenblk;
                xz.conservativeResize(lenxz,L);
                xz.block(lenxz_old, 0, xz.rows()-lenxz_old, L)= bl_BT_pad(n.segment(1, n.size()-1), f, x.block((i-1)*b_edge, 0, b_edge, L));
            }
        }
        return xz;
    }
}



VectorXcd BT_reconstruct(const VectorXi& n, const int& f, const VectorXcd& bz) // テスト済み
{
    //Reconstruction of b from bz
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     bz                   input vector
     /*************************** Output ***************************************************
     b                    output vector (trimmed)
     ***************************************************************************************/

    VectorXcd b;
    if(n.size() < 1){
        for(int i= f; i<=f*f ; i+= f){
            b.conservativeResize(b.size()+1);
            b(b.size()-1)= bz(i-1);
        }
        return b;
    }else{
        int b_edge= f*f*(2*n.segment(1, n.size()-1).array()-1).prod();
        int this_n= n(0);
        for(int i= this_n; i <= 2*this_n-1; ++i){
            int lenb_old= b.size();
            int lenblk= BT_reconstruct(n.segment(1, n.size()-1), f, bz.segment((i-1)*b_edge,b_edge)).size();
            int lenb= lenb_old+lenblk;
            b.conservativeResize(lenb);
            b.segment(lenb_old,lenb-lenb_old)= VectorXcd::Zero(lenb-lenb_old);
            b.segment(lenb_old,lenblk)= BT_reconstruct(n.segment(1, n.size()-1), f, bz.segment((i-1)*b_edge,b_edge));
        }
        return b;
    }
}


MatrixXcd bl_BT_reconstruct(const VectorXi& n, const int& f, const MatrixXcd& bz)
{
    //Reconstruction of b from bz
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     bz                   input (fft_length x L) matrix
     /*************************** Output ***************************************************
     b                    output (n x L) matrix (trimmed)
     ***************************************************************************************/

    MatrixXcd b;
    int L= bz.cols();
    if(n.size() < 1){
        for(int i= f; i<=f*f ; i+= f){
            b.conservativeResize(b.rows()+1,L);
            b.row(b.rows()-1)= bz.row(i-1);
        }
        return b;
    }else{
        int b_edge= f*f*(2*n.segment(1, n.size()-1).array()-1).prod();
        int this_n= n(0);
        for(int i= this_n; i <= 2*this_n-1; ++i){
            int lenb_old= b.rows();
            int lenblk= bl_BT_reconstruct(n.segment(1, n.size()-1), f, bz.block((i-1)*b_edge, 0, b_edge, L)).rows();
            int lenb= lenb_old+lenblk;
            b.conservativeResize(lenb,L);
            b.block(lenb_old,0,lenb-lenb_old,L).setZero();
            b.block(lenb_old,0,lenblk,L)= bl_BT_reconstruct(n.segment(1, n.size()-1), f, bz.block((i-1)*b_edge,0,b_edge,L));
        }
        return b;
    }
}


VectorXcd MBT_fft_init(const VectorXi& n, const int& f, const double& lf, const double& k0)
{
    //Prepare the fourier-transformed MBT projection of MBT matrix
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     lf                   physical length per lattice spacing (length scale factor)
     k0                   free space wavenumber
     /*************************** Output ***************************************************
     Au_til                result vector
     ***************************************************************************************/

    MatrixXi n_ind_init(2,n.size());
    n_ind_init= MatrixXi::Constant(2,n.size(),1);
    int fft_length= f*f*(2*n.array()-1).prod();

    int level_init= 1;
    //MBT projection of MBT matrix A;
    VectorXcd Au= BT_fft(n,f,n_ind_init,level_init,lf,k0);
    if(!(fft_length == Au.size())){
        cerr << "wrong Au size !" << endl;
        exit(EXIT_FAILURE);
    }

    Au=Au.reverse().eval();
    FFT<double> fft;
    VectorXcd Au_til= fft.fwd(Au);
    return Au_til;
}

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

MatrixXcd bl_MBT_fft_mvp(const VectorXi& n,const int& f, const VectorXcd& Au_til, const MatrixXcd& P_hat)
{
    //fast matrix-vector product using FFT. matrix is MBT
    /*************************** Input parameters *****************************************
     n=[n1,n2,...,nM]     number of BT_blocks at each level, where M is the number of BT levels.
     f                    size of final dense f*f matrix (e.g., f=3 for 3*3 matrix)
     Au_til               fourier-transformed MBT projection of MBT matrix A
     P_hat                input (n x L) matrix
     /*************************** Output ***************************************************
     Q_hat                result (n x L) matrix
     ***************************************************************************************/
     int fft_length= Au_til.size();
     int L= P_hat.cols();

     //MBT Projection of data vector
     MatrixXcd P_hatz= bl_BT_pad(n, f, P_hat);
     int len_bl_bt_pad= P_hatz.rows();
     P_hatz.conservativeResize(fft_length,L);
     P_hatz.block(len_bl_bt_pad,0,fft_length-len_bl_bt_pad,L).setZero();

     MatrixXcd P_hatz_til(fft_length,L);
     MatrixXcd Q_hatz_til(fft_length,L);
     MatrixXcd Q_hatz(fft_length,L);
     FFT<double> fft;
     for(int l= 0; l < L; ++l){
         P_hatz_til.col(l)= fft.fwd(P_hatz.col(l));
         Q_hatz_til.col(l)= Au_til.cwiseProduct(P_hatz_til.col(l));
         Q_hatz.col(l)= fft.inv(Q_hatz_til.col(l));
     }

     MatrixXcd Q_hat= bl_BT_reconstruct(n, f, Q_hatz);

     return Q_hat;
}



int main() {

    int f= 6; // ! f is the size of final dense f*f matrix (f=3 for CED, f=6 for CEMD)
    int nx= 10;
    int ny= 10;
    int nz= 15;
    VectorXi n(3);
    n << nx,ny,nz; // n=(nx,ny,nz), the dimension of cuboid domain
    cout << n << endl;



    double k0= 1;
    double lf= 1;


/*
    vector<Vector3d> pos_element;
    for(int ix= 0; ix < nx; ++ix){
        for(int iy= 0; iy < ny; ++iy){
            for(int iz= 0; iz < nz; ++iz){
                Vector3d pos_lattice(1.0*ix,1.0*iy,1.0*iz);
                pos_element.push_back(pos_lattice*lf);
            }
        }
    }

    MatrixXcd A(f*pos_element.size(),f*pos_element.size());
    //VectorXcd P(f*pos_element.size());

    cout << "matrix assingment ... " <<  endl;
    for(int m= 0; m < pos_element.size(); ++m){
        for(int k= 0; k < pos_element.size(); ++k){
            switch (f)
            {
            case 3:
                {
                    if(k == m){
                        A.block<3,3>(3*k,3*k)= Matrix3cd::Zero();
                    }else{
                        Vector3d Xkm= pos_element[k]-pos_element[m];
                        double Xkm_abs= Xkm.norm();
                        Matrix3d XkmXkm= Xkm*Xkm.transpose();
                        MatrixXcd Gkm;
                        Gkm= k0*k0*(MatrixXcd::Identity(3,3)-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
                        Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs/Xkm_abs)*(MatrixXcd::Identity(3,3)-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
                        Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs;
                        A.block<3,3>(3*k,3*m)= -Gkm;
                    }

                break;
                }
            case 6:
                {
                    if(k == m){
                        A.block<3,3>(6*k,6*k)= MatrixXcd::Zero(3,3);
                        A.block<3,3>(6*k+3,6*k+3)= MatrixXcd::Zero(3,3);
                    }else{
                        Vector3d Xkm= pos_element[k]-pos_element[m];
                        double Xkm_abs= Xkm.norm();
                        Vector3d ukm= Xkm/Xkm_abs;
                        Matrix3d ukm_vp= MatrixXd::Zero(3,3);
                        ukm_vp(0,1)= -ukm(2);
                        ukm_vp(0,2)= ukm(1);
                        ukm_vp(1,0)= ukm(2);
                        ukm_vp(1,2)= -ukm(0);
                        ukm_vp(2,0)= -ukm(1);
                        ukm_vp(2,1)= ukm(0);
                        Matrix3d XkmXkm= Xkm*Xkm.transpose();
                        Matrix3cd Gkm;
                        Gkm= k0*k0*(MatrixXcd::Identity(3,3)-complex<double>(1.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
                        Gkm= Gkm-((1.0-complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs/Xkm_abs)*(MatrixXcd::Identity(3,3)-complex<double>(3.0+0.0i)*XkmXkm/Xkm_abs/Xkm_abs);
                        Gkm= Gkm*exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs;
                        Matrix3cd Fkm;
                        Fkm= -k0*k0*(exp(complex<double>(0.0+1.0i)*k0*Xkm_abs)/Xkm_abs)*(1.0-1.0/(complex<double>(0.0+1.0i)*k0*Xkm_abs))*ukm_vp;
                        A.block<3,3>(6*k,6*m)= -Gkm;
                        A.block<3,3>(6*k,6*m+3)= -Fkm;
                        A.block<3,3>(6*k+3,6*m)= Fkm;
                        A.block<3,3>(6*k+3,6*m+3)= -Gkm;
                    }

                break;
                }
            default:

                break;
            }
        }
    }

    cout << "matrix assingment finished ... " <<  endl;
*/
    int fft_length= f*f*(2*n.array()-1).prod();

    FFT<double> fft;
    VectorXcd Au_til = MBT_fft_init(n, f, lf, k0);


    //length of data vector
    int vec_length= f*n.prod();


    //Vector before MVP x
    VectorXcd x(vec_length);
    for(int i=0; i< vec_length; ++i) x(i) = i+1;
    //cout << "vector before MVP x = " << endl;
    cout << "sum(x)= " << x.sum() << endl<< endl;

    VectorXcd b= MBT_fft_mvp(n,f, Au_til, x);

    cout << "size(b)= " << b.size() << endl;
    cout << "sum(b)= " << b.sum() << endl;



    int L=3;
    MatrixXcd X(vec_length,L);
    for(int i=0; i< vec_length; ++i) {
        for(int l=0; l<L; ++l){
            X(i,l) = (i+1);
        }
    }


    MatrixXcd B= bl_MBT_fft_mvp(n,f, Au_til, X);


    cout << "B.rows(), B.cols()= " << B.rows() << ", " << B.cols() << endl;
    for(int l= 0; l < L; ++l){
        cout << "B.col(l).sum= " << B.col(l).sum() << endl;
    }







}
