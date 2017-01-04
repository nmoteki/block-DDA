#ifndef INCLUDED_MATERIAL
#define INCLUDED_MATERIAL

#ifndef INCLUDED_MATERIAL_NAME
#include "Material_Name.hpp"
#endif

#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <sstream>

class Material { // material properties for target and surrounding medium

    std::string filename;
    /* format of material data
    ----------------------------
    wl_0[0], eper[0], mper[0]
    wl_0[1], eper[1], mper[1]
    ...    , ...    , ...
    ...    , ...    , ...
    ----------------------------
    format of eper and mper ==   (real_part,imaginary_part)
    */ //number of rows is arbitrary
    Material_Name material_name;
    double density; //density in gcm-3
    std::vector<double> wl0; // free_space_wavelength (um)
    std::vector<std::complex<double>> eper;
    std::vector<std::complex<double>> mper;

    std::complex<double> const_eper;
    std::complex<double> const_mper;
public:
    Material(){}; //default constructor
    Material(Material_Name material_name_set, double density_set, std::complex<double> const_eper_set, std::complex<double> const_mper_set); // constructor
    std::complex<double> get_eper(double free_space_wavelength) const; // return wavelength-interpolated eper (return const_eper if Material == CONSTANT)
    std::complex<double> get_mper(double free_space_wavelength) const; // return wavelength-interpolated mper (return const_mper if Material == CONSTANT)
    double get_density() const {return density;} // get density
    Material_Name get_name() const {return material_name;} // get material name
    void show_info() const;
};

#endif
