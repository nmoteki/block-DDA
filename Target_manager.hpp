#ifndef INCLUDED_TARGET_MANAGER
#define INCLUDED_TARGET_MANAGER

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <omp.h>

#include <Eigen/Dense>

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

#include "gmie.hpp"

class Target_manager {
    std::vector<Target> target_queue;
public:
    bool is_in_target_volume(const Eigen::Vector3d& pos_vec, Shape_Name shape_name, const std::vector<double>& side_lengths);
    void add_to_target_queue(Shape_Name shape_name, int nmax, const std::vector<double>& side_lengths, const Material& elem_material); // nmax:= max(nx,ny,nz)
    void add_to_target_queue(const std::string& target_file_name, const Material& material_0, const Material& material_1, double scale_factor); // read target data from file
    /* format of target file
    ---------------------------------------------
    aggregate_pp_radius
    coating_radius
    aggregate_volume
    coating_volume
    nx, ny, nz
    lf
    pos[0], vol[0], address[0], material_index[0]
    pos[1], vol[1], address[1], material_index[1]
    pos[2], vol[2], address[2], material_index[2]
    ...   , ...   , ...       , ...
    ...   , ...   , ...       , ...
    --------------------------------------------
    */

    int get_queue_size() const {return target_queue.size();};
    void set_optical_constant_at_this_wl0(Target& target, double free_space_wavelength);
    void set_polarizability(Target& target, int set_f,const Incident_Wave& incwave);
    void set_incident_fields_in_target(Target&, const Incident_Wave& incwave);

    Target pick_target(int target_index) {return target_queue[target_index];}

};

#endif
