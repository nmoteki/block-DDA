#ifndef INCLUDED_INCIDENT_WAVE
#include "Incident_Wave.hpp"
#endif

using namespace Eigen;
using namespace std;

Incident_Wave::Incident_Wave(const Material& surrounding_medium, const double free_space_wavelength) {
    wl0= free_space_wavelength;
    k0= 2*M_PI/wl0;
    eper_medium= surrounding_medium.get_eper(wl0);
    mper_medium= surrounding_medium.get_mper(wl0);
    k= k0*sqrt(eper_medium/mper_medium);
    uinc0 << 0.0,0.0,1.0; //default propagation vector (+z direction)
    einc0 << 1.0,0.0,0.0; //default E-polarization vector (+x direction)
    hinc0 << 0.0,sqrt(eper_medium/mper_medium),0.0; //corresponding H-polarization vector (+y direction)

    // default incident wave:
    // x-polarized E-field, propagaing into +z direction
    rot_mat.push_back(Matrix3d::Identity());

    uinc.push_back(uinc0);
    einc.push_back(einc0);
    hinc.push_back(hinc0);
    N_inc= 1;

}

void Incident_Wave::set_random_incident_directions(int L){
    // add Ninc-1 incident waves, random direction
    //default_random_engine re(random_device{}()); // random seed
    default_random_engine re; //default fixed seed
    uniform_real_distribution<double> urf {0,1.0};

    for(int i= 1; i < L; ++i){
        double theta= acos(2.0*urf(re)-1.0); // polar angle
        double phi= 2*M_PI*urf(re); // azimuth angle
        Vector3d nvec;
        nvec << sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta); //axis of rotation randomly choosen over unit sphere
        double psi= 2*M_PI*urf(re); // randomly choosen rotation angle
        Matrix3d mat;
        mat= AngleAxisd(psi,nvec);
        rot_mat.push_back(mat); // rotation matrix for i-th propagation direction
        uinc.push_back(rot_mat[i]*uinc0);
        einc.push_back(rot_mat[i]*einc0);
        hinc.push_back(rot_mat[i]*hinc0);
    }
    N_inc = L;
}

void Incident_Wave::show_info() const{
    cout << "------Incident_Wave Info------" << endl;
    cout << "free_space_wavelength == " << wl0 << endl;
    cout << "free_space_wavenumber == " << k0 << endl;
    cout << "eper_medium == " << eper_medium << endl;
    cout << "mper_medium == " << mper_medium << endl;
    cout << "in_medium_wavenumber == " << k << endl;
    cout << "N_inc == " << N_inc << endl;
    cout << "uinc == " << endl;
    for(const auto& x : uinc)  cout << x.transpose() << endl;
    cout << "------------------------------" << endl << endl;
}
