#ifndef INCLUDED_INCIDENT_WAVE
#define INCLUDED_INCIDENT_WAVE

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <random>
#include <Eigen/Dense>

#ifndef INCLUDED_MATERIAL
#include "Material.hpp"
#endif


class Incident_Wave {
    Eigen::Vector3d uinc0; //default propagation vector (+z direction)
    Eigen::Vector3cd einc0; //x-polarized E amplitude
    Eigen::Vector3cd hinc0; //H amplitude corresponding to einc0
    double wl0; // free-space wavelength (=c/w)
    double k0; // free-space wavenumber
    std::complex<double> eper_medium;
    std::complex<double> mper_medium;
    std::complex<double> k; // wavenumber in surrounding medium
    int N_inc;
    std::vector<Eigen::Matrix3d> rot_mat; // std::vector of 3D rotation matrices for individual propagation directions
    std::vector<Eigen::Vector3d> uinc; // unit vector of propagation directions
    std::vector<Eigen::Vector3cd> einc; // unit vector of propagation directions
    std::vector<Eigen::Vector3cd> hinc; // unit vector of propagation directions
public:
    Incident_Wave(){}; //default constructor
    Incident_Wave(const Material& surrounding_medium, const double free_space_wavelength);
    double get_wl0() const {return wl0;}
    std::complex<double> get_eper() const {return eper_medium;}
    std::complex<double> get_mper() const {return mper_medium;}
    std::complex<double> get_k() const {return k;}
    int get_N_inc() const {return N_inc;}
    Eigen::Vector3d get_uinc(int l) const {return uinc[l];}
    Eigen::Vector3cd get_einc(int l) const {return einc[l];}
    Eigen::Vector3cd get_hinc(int l) const {return hinc[l];}
    void set_random_incident_directions(int L);
    void show_info() const;
};

#endif
