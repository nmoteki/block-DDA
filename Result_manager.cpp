#ifndef INCLUDED_RESULT_MANAGER
#include "Result_manager.hpp"
#endif


using namespace std;
using namespace Eigen;


void Result_manager::add_to_result_queue(const Target& target, const Incident_Wave& incwave, const Material& medium, const int wl0_index, const int target_index){
    double um3_to_cm3 {1.0e-12};

    Result result;
    result.wl0_index= wl0_index;
    result.target_index= target_index;
    result.wl0= incwave.get_wl0();
    result.N_inc= incwave.get_N_inc();
    result.material_name_medium= medium.get_name();
    result.eper_medium= medium.get_eper(result.wl0);
    result.mper_medium= medium.get_mper(result.wl0);

    result.target_fname= target.target_fname;
    auto it = find_if(target.eper.begin(), target.eper.end(), [=](complex<double> eper){ return (eper != target.eper[0]);});
    auto index_begin_second_material= distance(target.eper.begin(),it);
    if(it==target.eper.end()){ // homogeneous material
        result.is_homogeneous= true;
        result.material_name_particle= target.material[0].get_name();
        result.eper_particle= target.eper[0];
        result.mper_particle= target.mper[0];
        result.volume_particle= accumulate(target.vol.begin(),target.vol.end(),0.0);
        result.density_particle= target.material[0].get_density();
        result.mass_particle= result.volume_particle*um3_to_cm3*result.density_particle;
        result.ve_radius_particle= cbrt(result.volume_particle*3.0/(4.0*M_PI));
        result.ve_cross_section_particle= M_PI*result.ve_radius_particle*result.ve_radius_particle;
    }else{ //internal mixture of two materials
        result.is_homogeneous= false;
        result.material_name_core= target.material[0].get_name();
        result.eper_core= target.eper[0];
        result.mper_core= target.mper[0];
        result.volume_core= accumulate(target.vol.begin(),target.vol.begin()+index_begin_second_material,0.0);
        result.density_core= target.material[0].get_density();
        result.mass_core= result.volume_core*um3_to_cm3*result.density_core;
        result.ve_radius_core= cbrt(result.volume_core*3.0/(4.0*M_PI));
        result.ve_cross_section_core= M_PI*result.ve_radius_core*result.ve_radius_core;

        result.material_name_coat= target.material[index_begin_second_material].get_name();
        result.eper_coat= target.eper[index_begin_second_material];
        result.mper_coat= target.mper[index_begin_second_material];
        result.volume_coat= accumulate(target.vol.begin()+index_begin_second_material,target.vol.end(),0.0);
        result.coat_to_core_vratio= result.volume_coat/result.volume_core;
    }


    int L= target.EHinc.cols();
    result.Cext.resize(L);
    result.Cabs.resize(L);
    result.Csca.resize(L);
    result.Qext.resize(L);
    result.Qabs.resize(L);
    result.Qsca.resize(L);
    result.SSA.resize(L);
    result.MAC.resize(L);


    complex<double> k= incwave.get_k();
    switch (target.f)
    {
    case 3: // electric target
        {
            for(int l= 0; l < L; ++l){
                if(imag(k)==0.0){
                    result.Cext[l]= 4*M_PI*real(k)*imag(target.EHinc.col(l).dot(target.PM.col(l)));
                    result.Cabs[l]= 4*M_PI*real(k)*(imag(conj(1.0/target.alpha_E(0)))-2.0*real(k)*real(k)*real(k)/3.0)*target.PM.col(l).squaredNorm();
                }else{
                    cerr << "non zero imag(k) is not supported" << endl;
                    exit(EXIT_FAILURE);
                }
                result.Qext[l]= result.Cext[l]/result.ve_cross_section_particle;
                result.Qabs[l]= result.Cabs[l]/result.ve_cross_section_particle;
                result.Csca[l]= result.Cext[l]-result.Cabs[l];
                result.Qsca[l]= result.Qext[l]-result.Qabs[l];
                result.SSA[l]= (result.Cext[l]-result.Cabs[l])/result.Cext[l];

                double um2_to_m2 {1.0e-12};
                double mass= result.is_homogeneous ? result.mass_particle : result.mass_core;
                result.MAC[l]= result.Cabs[l]*um2_to_m2/mass;
            }
            break;
        }
    case 6: // magnetic target
        {
            for(int l= 0; l < L; ++l){
                if(imag(k)==0.0){
                    result.Cext[l]= 4*M_PI*real(k)*imag(target.EHinc.col(l).dot(target.PM.col(l)));
                    result.Cabs[l]= 0.0;
                    for(int j= 0; j < target.num_element_occupy; ++j){
                        result.Cabs[l]+= (imag(conj(1.0/target.alpha_E(0)))-2.0*real(k)*real(k)*real(k)/3.0)*target.PM.block<3,1>(6*j,l).squaredNorm(); // contribution of P
                        result.Cabs[l]+= (imag(conj(1.0/target.alpha_H(0)))-2.0*real(k)*real(k)*real(k)/3.0)*target.PM.block<3,1>(6*j+3,l).squaredNorm(); // contribution of M
                    }
                    result.Cabs[l]*= 4*M_PI*real(k);
                }else{
                    cerr << "non zero imag(k) is not supported" << endl;
                    exit(EXIT_FAILURE);
                }

                result.Qext[l]= result.Cext[l]/result.ve_cross_section_particle;
                result.Qabs[l]= result.Cabs[l]/result.ve_cross_section_particle;
                result.Csca[l]= result.Cext[l]-result.Cabs[l];
                result.Qsca[l]= result.Qext[l]-result.Qabs[l];
                result.SSA[l]= (result.Cext[l]-result.Cabs[l])/result.Cext[l];

                double um2_to_m2 {1.0e-12};
                double mass= result.is_homogeneous ? result.mass_particle : result.mass_core;
                result.MAC[l]= result.Cabs[l]*um2_to_m2/mass;
            }
            break;
        }
    default:
        cerr << "unknown f value" << endl;
        exit(EXIT_FAILURE);
        break;
    }

    result.f= target.f;
    result.nx= target.n(0);
    result.ny= target.n(1);
    result.nz= target.n(2);
    result.num_element_occupy= target.num_element_occupy;
    result.num_element_cuboid= target.n.prod();
    result.DDA_solver_time= target.DDA_solver_time;

    result_queue.push_back(result);
}



tuple<double,double,double,double> Result_manager::volume_equivalent_Mie_result(const Result& result){
    // Exact solution for sphere
    double um2_to_m2 {1.0e-12};
    double wl0= result.wl0;
    complex<double> eper_m= result.eper_medium;
    complex<double> mper_m= result.mper_medium;
    complex<double> eper_p;
    complex<double> mper_p;
    double mass;
    double r_p;
    double geom_cs;
    if(result.is_homogeneous){
        eper_p= result.eper_particle;
        mper_p= result.mper_particle;
        mass= result.mass_particle;
        r_p= result.ve_radius_particle;
        geom_cs= result.ve_cross_section_particle;
    }else{
        eper_p= result.eper_core;
        mper_p= result.mper_core;
        mass= result.mass_core;
        r_p= result.ve_radius_core;
        geom_cs= result.ve_cross_section_core;
    }

    int nang {3}; // number of scattering angles between 0-180 deg for scattering matrix element S1 S2
    double Qsca_mie, Qabs_mie, Qext_mie;
    vector<complex<double>> S1,S2;
    tie(Qsca_mie,Qabs_mie,Qext_mie,S1,S2)= gmie(wl0, r_p, nang, eper_p, mper_p, eper_m, mper_m);

    double Cext_mie= Qext_mie*M_PI*r_p*r_p;
    double Cabs_mie= Qabs_mie*M_PI*r_p*r_p;
    double Csca_mie= Qsca_mie*M_PI*r_p*r_p;

    double MAC_mie= Cabs_mie*um2_to_m2/mass;

    return forward_as_tuple(Cext_mie,Cabs_mie,Csca_mie,MAC_mie);

}

//define << for output single result
ostream& operator<<(ostream& os, const Result& result){

    Basic_Stats<double> Csca_stats(result.Csca);
    Basic_Stats<double> MAC_stats(result.MAC);
    Basic_Stats<double> SSA_stats(result.SSA);

    os << "####################### Start of result section #######################" << endl;
    os << "wl0_index" << '\t' << result.wl0_index << endl;
    os << "target_index" << '\t' << result.target_index << endl;
    os << endl;
    os << "----------------- Incident wave properties -------------------" << endl;
    os << "wl0" << '\t' << result.wl0 << endl
       << "N_inc" << '\t' << result.N_inc << endl;
    os << endl;
    os << "-------------------- Medium properties -----------------------" << endl;
    os << "eper_medium"<< '\t' << result.eper_medium << endl
       << "mper_medium"<< '\t' << result.mper_medium << endl;
    os << endl;

    os << "-------------------- Target properties ----------------------" << endl;
    os << "target_fname" << '\t' << result.target_fname << endl;
    os << "is_homogeneous" << '\t' << result.is_homogeneous << endl;
    if(result.is_homogeneous){
        os //<< "material_name_particle" << '\t' << result.material_name_particle << endl
           << "eper_particle"<< '\t' << result.eper_particle << endl
           << "mper_particle"<< '\t' << result.mper_particle << endl
           << "volume_particle"<< '\t' << result.volume_particle << endl
           << "mass_particle"<< '\t' << result.mass_particle << endl
           << "density_particle"<< '\t' << result.density_particle << endl
           << "ve_radius_particle" << '\t' << result.ve_radius_particle << endl;
    }else{
        os //<< "material_name_core" << '\t' << result.material_name_core << endl
           << "eper_core"<< '\t' << result.eper_core << endl
           << "mper_core"<< '\t' << result.mper_core << endl
           << "volume_core"<< '\t' << result.volume_core << endl
           << "mass_core"<< '\t' << result.mass_core << endl
           << "density_core"<< '\t' << result.density_core << endl
           << "ve_radius_core" << '\t' << result.ve_radius_core << endl
           //<< "material_name_coat" << '\t' << result.material_name_coat << endl
           << "eper_coat"<< '\t' << result.eper_coat << endl
           << "mper_coat"<< '\t' << result.mper_coat << endl
           << "volume_coat"<< '\t' << result.volume_coat << endl
           << "coat_to_core_vratio" << '\t' << result.coat_to_core_vratio << endl;

    }
    os << endl;
    os << "----------------- Computed optical properties -----------------" << endl;
    os << "mean_Csca" << '\t' << Csca_stats.mean() << endl
       << "std_Csca" << '\t' << Csca_stats.std() << endl
       << "mean_MAC" << '\t' << MAC_stats.mean() << endl
       << "std_MAC" << '\t' << MAC_stats.std() << endl
       << "mean_SSA" << '\t' << SSA_stats.mean() << endl
       << "std_SSA" << '\t' << SSA_stats.std() << endl;

    os << endl;
       os << "------------------- Computational info --------------------" << endl;
       os << "f" << '\t' << result.f << endl;
       os << "nx" << '\t' << result.nx << endl;
       os << "ny" << '\t' << result.ny << endl;
       os << "nz" << '\t' << result.nz << endl;
       os << "num_element_cuboid" << '\t' << result.num_element_cuboid << endl;
       os << "num_element_occupy" << '\t' << result.num_element_occupy << endl;
       os << "DDA_solver_time" << '\t' << result.DDA_solver_time << endl
          << "DDA_solver_time_per_Ninc" << '\t' << result.DDA_solver_time/result.N_inc << endl;
    os << endl;
    os << "************************ End of result section ***********************" << endl;

}
