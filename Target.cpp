#ifndef INCLUDED_TARGET
#include "Target.hpp"
#endif

using namespace std;

double Target::vol_eq_radius() const{
    double total_vol= accumulate(vol.begin(),vol.end(),0.0);
    double ve_radius= cbrt(total_vol*3.0/(4.0*M_PI));
    return ve_radius;
}

double Target::vol_eq_radius_element(int index) const{
    if(index < vol.size()) {
        return  cbrt(vol[index]*3.0/(4.0*M_PI));
    }
    cerr << " incorrect index " << endl; // To do "out_of_range"
    exit(EXIT_FAILURE);
}

void Target::show_info() const {
    cout << "-------- Target info --------" << endl;
    cout << "fname== " << target_fname << endl;
    cout << "n== " << n.transpose() << endl;
    cout << "f== " << f << endl;
    cout << "num_element_occupy== " << num_element_occupy << endl;
    cout << "eper[0]== " << eper[0] << endl;
    cout << "mper[0]== " << mper[0] << endl;
    cout << "-----------------------------" << endl << endl;
}
