#ifndef INCLUDED_RESULT_MANAGER
#define INCLUDED_RESULT_MANAGER

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ostream>
#include <tuple>
#include <Eigen/Dense>

#ifndef INCLUDED_INCIDENT_WAVE
#include "Incident_Wave.hpp"
#endif
#ifndef INCLUDED_TARGET
#include "Target.hpp"
#endif
#ifndef INCLUDED_TARGET_MANAGER
#include "Target_manager.hpp"
#endif
#ifndef INCLUDED_BASIC_STATS
#include "Basic_Stats.hpp"
#endif
#ifndef INCLUDED_MATERIAL_NAME
#include "Material_Name.hpp"
#endif

#include "gmie.hpp"


struct Result {
    int wl0_index;
    int target_index;

    double wl0;
    double N_inc;

    Material_Name material_name_medium;
    complex<double> eper_medium;
    complex<double> mper_medium;

    std::string target_fname;

    bool is_homogeneous; // true if target is homogeneous material

    Material_Name material_name_particle;
    complex<double> eper_particle;
    complex<double> mper_particle;
    double volume_particle; //um3
    double mass_particle;
    double density_particle;
    double ve_radius_particle; //um
    double ve_cross_section_particle; //um2

    Material_Name material_name_core;
    complex<double> eper_core;
    complex<double> mper_core;
    double volume_core; //um3
    double mass_core; //g
    double density_core;
    double ve_radius_core; //um
    double ve_cross_section_core; //um2

    Material_Name material_name_coat;
    complex<double> eper_coat;
    complex<double> mper_coat;
    double volume_coat; //um3
    double coat_to_core_vratio;

    std::vector<double> Cext;
    std::vector<double> Cabs;
    std::vector<double> Csca;
    std::vector<double> Qext;
    std::vector<double> Qabs;
    std::vector<double> Qsca;
    std::vector<double> SSA; // Csca/Cext
    std::vector<double> MAC; // mass_absorption_cross_section of core (m2/g)

    int f;
    int nx;
    int ny;
    int nz;
    int num_element_cuboid;
    int num_element_occupy;
    double DDA_solver_time;
};

class Result_manager {
    std::vector<Result> result_queue;
public:
    void add_to_result_queue(const Target& target, const Incident_Wave& incwave, const Material& medium, const int wl0_index, const int target_index);
    Result pick_result(int result_index) {return result_queue[result_index];}
    Result pick_result(const Target_manager& target_manager, int wl0_index, int target_index){
        int result_index= wl0_index*target_manager.get_queue_size()+target_index;
        return result_queue[result_index];
    }
    //void compute_total_cross_sections(const Target&, const Incident_Wave&);
    tuple<double,double,double,double> volume_equivalent_Mie_result(const Result& result);
};

std::ostream& operator<<(std::ostream& os, const Result& result);

#endif
