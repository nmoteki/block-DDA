#ifndef INCLUDED_TARGET_MANAGER
#include "Target_manager.hpp"
#endif

using namespace std;
using namespace Eigen;

bool Target_manager::is_in_target_volume(const Vector3d& pos_vec, Shape_Name shape_name, const vector<double>& side_lengths){
    switch (shape_name)
    {
    case Shape_Name::SPHERE:
        {
            double r= side_lengths[0]/2.0; // sphere radius
            Vector3d centroid_pos_vec(r,r,r); // centroid position vector
            double dist_from_centroid= (pos_vec-centroid_pos_vec).norm();
            bool is_in_target_volume= (dist_from_centroid <= r);
            return is_in_target_volume;
            break;
        }
    case Shape_Name::CUBE: case Shape_Name::CUBOID:
        {
            bool is_in_target_volume= (pos_vec(0)<=side_lengths[0] && pos_vec(1)<=side_lengths[1] && pos_vec(2)<=side_lengths[2]);
            return is_in_target_volume;
            break;
        }
    case Shape_Name::ELLIPSOID:
        {
            double a= side_lengths[0]/2.0; // radius along x axis
            double b= side_lengths[1]/2.0; // radius along y axis
            double c= side_lengths[2]/2.0; // radius along z axis
            Vector3d centroid_pos_vec(a,b,c); // centroid position vector
            double ellip_func= (pos_vec(0)-centroid_pos_vec(0))*(pos_vec(0)-centroid_pos_vec(0))/(a*a)
                                +(pos_vec(1)-centroid_pos_vec(1))*(pos_vec(1)-centroid_pos_vec(1))/(b*b)
                                +(pos_vec(2)-centroid_pos_vec(2))*(pos_vec(2)-centroid_pos_vec(2))/(c*c);
            bool is_in_target_volume= (ellip_func<=1.0);
            return is_in_target_volume;
            break;
        }
    case Shape_Name::CYLINDER:
        {
            //for incident_wave propagating into +z direction, head_on orientation is assumed
            double r= side_lengths[0]/2.0; // cylinder radius
            double l= side_lengths[2]; // cylinder length

            double dist_from_axis= sqrt((pos_vec(0)-r)*(pos_vec(0)-r)+(pos_vec(1)-r)*(pos_vec(1)-r));
            bool is_in_target_volume= (dist_from_axis <= r && pos_vec(2) <=l);
            return is_in_target_volume;
            break;
        }
    case Shape_Name::ROUNDED_CYLINDER:
        {
            //for incident_wave propagating into +z direction, head_on orientation is assumed
            double r= side_lengths[0]/2.0; // cylinder radius
            double l= side_lengths[2]; // cylinder length (top-to-top length)

            double dist_from_axis= sqrt((pos_vec(0)-r)*(pos_vec(0)-r)+(pos_vec(1)-r)*(pos_vec(1)-r)); //(flat part)
            double dist_from_center_z_ls_r= sqrt((pos_vec(0)-r)*(pos_vec(0)-r)+(pos_vec(1)-r)*(pos_vec(1)-r)+(pos_vec(2)-r)*(pos_vec(2)-r));
            double dist_from_center_z_gt_lplusr= sqrt((pos_vec(0)-r)*(pos_vec(0)-r)+(pos_vec(1)-r)*(pos_vec(1)-r)+(pos_vec(2)-(l-r))*(pos_vec(2)-(l-r)));
            bool is_in_target_volume;

            if(pos_vec(2) <= r){
                is_in_target_volume= (dist_from_center_z_ls_r <= r);
                return is_in_target_volume;
            }
            if(r < pos_vec(2) || pos_vec(2) < l-r){
                is_in_target_volume= (dist_from_axis <= r);
                return is_in_target_volume;
            }
            if(pos_vec(2) >= l-r){
                is_in_target_volume= (dist_from_center_z_gt_lplusr <= r);
                return is_in_target_volume;
            }
            break;
        }
    default:
        {
            cerr << "unknown Shape_Name " << endl;
            exit(EXIT_FAILURE);
            break;
        }
    }
}

void Target_manager::add_to_target_queue(Shape_Name shape_name, int nmax, const vector<double>& side_lengths, const Material& elem_material){
    Target target; // new target
    if(!(side_lengths.size()==3)) {
        cerr << "wrong length of vector<double> side_lengths" << endl;
        exit(EXIT_FAILURE);
    }
    double max_side_length= *max_element(side_lengths.begin(),side_lengths.end());
    target.lf= max_side_length/nmax; // initialize target.lf
    target.n.resize(3);
    for(int index=0; index<3; ++index) target.n(index)= ceil(nmax*(side_lengths[index]/max_side_length)); // initialize target.n

    for(int ix= 0; ix < target.n(0); ++ix){
        for(int iy= 0; iy < target.n(1); ++iy){
            for(int iz= 0; iz < target.n(2); ++iz){
                int lattice_address= target.n(1)*target.n(2)*ix+target.n(2)*iy+iz;
                Vector3d pos_vec(target.lf*ix,target.lf*iy,target.lf*iz);
                if(is_in_target_volume(pos_vec,shape_name,side_lengths)){
                    target.pos.push_back(pos_vec);
                    target.address.push_back(lattice_address);
                    double elem_vol= target.lf*target.lf*target.lf;
                    target.vol.push_back(elem_vol);
                    target.material.push_back(elem_material);
                }
            }
        }
    }
    target.num_element_occupy= target.address.size();
    target_queue.push_back(target);
}

void Target_manager::add_to_target_queue(const string& target_file_name, const Material& material_0, const Material& material_1, double scale_factor){
    /* format of target file
    ---------------------------------------------
    aggregate_pp_radius
    coating_radius
    aggregate_volume
    coating_volume
    nx, ny, nz
    lf
    address[0], vol[0], material_index[0]
    address[1], vol[1], material_index[1]
    ...   , ...   , ...       , ...
    ...   , ...   , ...       , ...
    --------------------------------------------
    */

    Target target; // new target

    ifstream reading_file(target_file_name, ios::in);
    if(!reading_file){
        cerr << "file " << target_file_name << "cannot open" << endl;
        exit(EXIT_FAILURE);
    }
    cout << "reading " << target_file_name << "..." << endl;
    string line;
    double elem_vol;
    int lattice_address;
    int material_index;
    getline(reading_file,line);
    double pp_radius= stod(line);
    getline(reading_file,line);
    double aggregate_volume= stod(line);
    getline(reading_file,line);
    double coating_volume= stod(line);
    getline(reading_file,line);
    istringstream iss(line);
    target.n.resize(3);
    if(!(iss >> target.n(0) >> target.n(1) >> target.n(2))){
        cerr << "wrong data format in " << target_file_name << endl;
        exit(EXIT_FAILURE);
    }

    getline(reading_file,line);
    target.lf= stod(line)*scale_factor;
    //cout << "target.lf== " << target.lf << endl;

    while(getline(reading_file,line)){
        istringstream iss(line);

        if(!(iss >> lattice_address >> elem_vol  >> material_index)) {
            cerr << "wrong data format in " << target_file_name << endl;
            exit(EXIT_FAILURE);
        }
        int ix= lattice_address/(target.n(1)*target.n(2));
        int iy= (lattice_address-(target.n(1)*target.n(2))*ix)/target.n(2);
        int iz= lattice_address-(target.n(1)*target.n(2))*ix-target.n(2)*iy;

        Vector3d xyz(target.lf*ix,target.lf*iy,target.lf*iz);
        elem_vol *= (scale_factor*scale_factor*scale_factor);
        //cout << "elem_vol== " << elem_vol << endl;
        //cout << "lattice_address== " << lattice_address << endl;
        target.pos.push_back(xyz);
        target.vol.push_back(elem_vol);
        target.address.push_back(lattice_address);
        switch (material_index)
        {
        case 0:
            {
                target.material.push_back(material_0);
                break;
            }
        case 1:
            {
                target.material.push_back(material_1);
                break;
            }
        default:
            {
                cerr << "unknown material " << endl;
                exit(EXIT_FAILURE);
                break;
            }
        }
    }
    reading_file.close();
    target.target_fname= target_file_name;
    target.num_element_occupy= target.address.size();
    target_queue.push_back(target);
}

void Target_manager::set_optical_constant_at_this_wl0(Target& target, double free_space_wavelength){
    target.eper.clear();
    target.mper.clear();
    for(int i= 0; i<target.num_element_occupy; ++i){
        complex<double> eper= target.material[i].get_eper(free_space_wavelength);
        complex<double> mper= target.material[i].get_mper(free_space_wavelength);
        target.eper.push_back(eper);
        target.mper.push_back(mper);
    }
}

void Target_manager::set_polarizability(Target& target, int set_f, const Incident_Wave& incwave){
    target.f= set_f;
    double wl0= incwave.get_wl0();
    complex<double> eper_m= incwave.get_eper();
    complex<double> mper_m= incwave.get_mper();
    complex<double> k= incwave.get_k();
    switch (target.f)
    {
    case 3:
        {
            target.alpha_E.resize(target.num_element_occupy);
            for(int index= 0; index< target.num_element_occupy; ++index){
                double r_p= target.vol_eq_radius_element(index);
                complex<double> eper_p= target.eper[index];
                complex<double> mper_p= target.mper[index];
                vector<complex<double>> a,b,c,d;
                tie(a,b,c,d)= gmie_coeff(wl0, r_p, eper_p, mper_p, eper_m, mper_m);
                target.alpha_E(index)= 3.0*complex<double>(0.0+1.0i)*a[0]/(2.0*k*k*k);
            }
            break;
        }
    case 6:
        {
            target.alpha_E.resize(target.num_element_occupy);
            target.alpha_H.resize(target.num_element_occupy);
            for(int index= 0; index< target.num_element_occupy; ++index){
                double r_p= target.vol_eq_radius_element(index);
                complex<double> eper_p= target.eper[index];
                complex<double> mper_p= target.mper[index];
                vector<complex<double>> a,b,c,d;
                tie(a,b,c,d)= gmie_coeff(wl0, r_p, eper_p, mper_p, eper_m, mper_m);
                target.alpha_E(index)= 3.0*complex<double>(0.0+1.0i)*a[0]/(2.0*k*k*k);
                target.alpha_H(index)= 3.0*complex<double>(0.0+1.0i)*b[0]/(2.0*k*k*k);
            }

            break;
        }
    default:
        {
            cerr << "unknown target.f value " << endl;
            exit(EXIT_FAILURE);
            break;
        }
    }

}

void Target_manager::set_incident_fields_in_target(Target& target, const Incident_Wave& incwave) {
    int L= incwave.get_N_inc();
    complex<double> k= incwave.get_k();
    switch (target.f)
    {
    case 3:
        {
            target.EHinc.resize(3*target.num_element_occupy,L);
            for(int m= 0; m < target.num_element_occupy; ++m){
                for(int l= 0; l < L; ++l){
                    Vector3d uinc= incwave.get_uinc(l);
                    complex<double> phase_inc= exp(complex<double>(0.0+1.0i)*k*uinc.dot(target.pos[m]));
                    Vector3cd einc= incwave.get_einc(l);
                    target.EHinc.block<3,1>(3*m,l)= einc*phase_inc;
                }
            }
            break;
        }
    case 6:
        {
            target.EHinc.resize(6*target.num_element_occupy,L);
            for(int m= 0; m < target.num_element_occupy; ++m){
                for(int l= 0; l < L; ++l){
                    Vector3d uinc= incwave.get_uinc(l);
                    complex<double> phase_inc= exp(complex<double>(0.0+1.0i)*k*uinc.dot(target.pos[m]));
                    Vector3cd einc= incwave.get_einc(l);
                    Vector3cd hinc= incwave.get_hinc(l);
                    target.EHinc.block<6,1>(6*m,l) << einc*phase_inc, hinc*phase_inc;
                }
            }
            break;
        }
    default:
        cerr << "unknown target.f value " << endl;
        exit(EXIT_FAILURE);
        break;
    }
}
