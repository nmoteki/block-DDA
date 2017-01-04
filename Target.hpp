#ifndef INCLUDED_TARGET
#define INCLUDED_TARGET

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <complex>
#include <Eigen/Dense>

#ifndef INCLUDED_MATERIAL
#include "Material.hpp"
#endif


struct Target {
    Eigen::VectorXi n; // lattice size  n(0):=nx, n(1):=ny, n(2):=nz
    double lf; // length_scale factor := physical length of lattice spacing
    std::vector<Eigen::Vector3d> pos; // physical coordinate (x,y,z) of element = (ix,iy,iz)*lf
    std::vector<double> vol; // physical volume of element
    std::vector<int> address; // lattice address of element := ny*nz*(ix-1)+nz*(iy-1)+iz;
    std::vector<Material> material; // material of element
    std::vector<std::complex<double>> eper; // electric permittivity
    std::vector<std::complex<double>> mper; // magnetic permiability
    int num_element_occupy; // number of elements

    int f; // final dense MBT block size
    Eigen::VectorXcd alpha_E; // electric polarizability of element
    Eigen::VectorXcd alpha_H; // electric polarizability of element
    Eigen::MatrixXcd EHinc; // incident electric (and magnetic) fields in target
    Eigen::MatrixXcd PM; // electric (and magnetic) polarization

    std::string target_fname;

    double DDA_solver_time;

    double vol_eq_radius() const;
    double vol_eq_radius_element(int index) const;
    void show_info() const;
};

#endif
