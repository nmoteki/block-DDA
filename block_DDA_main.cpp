#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <tuple>
#include <string>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "gmie.hpp"

#ifndef INCLUDED_MATERIAL_NAME
#include "Material_Name.hpp"
#endif
#ifndef INCLUDED_SHAPE_NAME
#include "Shape_Name.hpp"
#endif
#ifndef INCLUDED_MATERIAL
#include "Material.hpp"
#endif
#ifndef INCLUDED_INCIDENT_WAVE
#include "Incident_Wave.hpp"
#endif
#ifndef INCLUDED_TARGET
#include "Target.hpp"
#endif
#ifndef INCLUDED_TARGET_MANAGER
#include "Target_manager.hpp"
#endif
#ifndef INCLUDED_SOLVER
#include "Solver.hpp"
#endif
#ifndef INCLUDED_RESULT_MANAGER
#include "Result_manager.hpp"
#endif
#ifndef INCLUDED_BASIC_STATS
#include "Basic_Stats.hpp"
#endif

using namespace Eigen;
using namespace std;



int main(int argc, char* argv[]){


int set_f{6}; // CED(3) or CEMD(6)

// prepare incident wave
double density_set= 1.0;
complex<double> const_eper_set= {1.0,0.0};
complex<double> const_mper_set= {1.0,0.0};
Material medium_material(Material_Name::VACUUM, density_set, const_eper_set, const_mper_set);
medium_material.show_info();


int N_inc{4}; // number of incident directions
//vector<double> wl0_list{0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5};
vector<double> wl0_list{0.3};
const_eper_set= {2.25,0.0};
Material coating_material(Material_Name::CONSTANT, density_set, const_eper_set, const_mper_set);


double core_density= {1.8}; //BC
Material core_material(Material_Name::BC, core_density, const_eper_set, const_mper_set);
double scale_factor{20.0e-3}; // pp_radius in um (BC)

//double core_density= {5.17}; //MAGNETITE
//double scale_factor{40.0e-3}; // pp_radius in um (MAGNETITE)
//Material core_material(Material_Name::MAGNETITE, core_density, const_eper_set, const_mper_set); // MAGNETITE


// prepare target
Target_manager target_manager;

string target_file_name;

//target_file_name= {"DDA_target_2_agg_N16384_kf1.0_Df2.8.out"}; // N=16384 n== (91, 88, 113) が MacPro 64GB における上限 一回 4.5hours
target_file_name= {"DDA_target_0_agg_N8_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N16_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N32_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N64_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N128_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N256_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
/*
target_file_name= {"DDA_target_0_agg_N512_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N1024_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N2048_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N4096_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N8192_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_0_agg_N16384_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N8_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N16_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N32_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N64_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N128_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N256_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N512_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N1024_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N2048_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N4096_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N8192_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);
target_file_name= {"DDA_target_2_agg_N16384_kf1.0_Df2.8.out"};
target_manager.add_to_target_queue(target_file_name, core_material, coating_material, scale_factor);

*/

Result_manager result_manager;

for(int wl0_index=0; wl0_index < wl0_list.size(); ++wl0_index){
    cout << endl << endl;
    double wl0_set= wl0_list[wl0_index];
    cout << "wl0_set== " << wl0_set << endl;
    Incident_Wave incident_wave(medium_material, wl0_set);
    incident_wave.set_random_incident_directions(N_inc);

    for(int target_index=0; target_index< target_manager.get_queue_size(); ++target_index){
        cout << endl << endl;
        cout << "target index== " << target_index << endl;
        Target target= target_manager.pick_target(target_index);
        target_manager.set_optical_constant_at_this_wl0(target, wl0_set);
        target_manager.set_polarizability(target, set_f, incident_wave);
        target_manager.set_incident_fields_in_target(target, incident_wave);
        target.show_info();

        // DDA solver
        Solver solver;
        solver.prepare_DDA_matrix(target, incident_wave);
        solver.show_info();
        solver.DDA_solve(target);

        result_manager.add_to_result_queue(target,incident_wave,medium_material,wl0_index,target_index);

    }
}

// result output
string ofilename{"DDA_result.dat"};
ofstream writing_file(ofilename, ios::out);
//ofstream writing_file(ofilename, ios::app);
for(int wl0_index=0; wl0_index < wl0_list.size(); ++wl0_index){
    for(int target_index=0; target_index< target_manager.get_queue_size(); ++target_index){

        Result result= result_manager.pick_result(target_manager,wl0_index,target_index);

        double Cext_mie,Cabs_mie,Csca_mie,MAC_mie;
        tie(Cext_mie,Cabs_mie,Csca_mie,MAC_mie)= result_manager.volume_equivalent_Mie_result(result);

        cout << result;
        cout << "MAC_mie" << '\t' << MAC_mie << endl;

        writing_file  << result;
        writing_file << "MAC_mie" << '\t' << MAC_mie << endl << endl;

    }
}


}
