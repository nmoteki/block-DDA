#include <Eigen/Dense>
using namespace Eigen;
using namespace std;


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
