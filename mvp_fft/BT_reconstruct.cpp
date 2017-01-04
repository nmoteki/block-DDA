#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

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
