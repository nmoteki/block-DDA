#ifndef INCLUDED_MATERIAL
#include "Material.hpp"
#endif

#ifndef INCLUDED_MATERIAL_NAME
#include "Material_Name.hpp"
#endif

using namespace std;

Material::Material(Material_Name material_name_set, double density_set, complex<double> const_eper_set={1.0,0.0}, complex<double> const_mper_set={1.0,0.0}){
    switch (material_name_set)
    {
    case Material_Name::VACUUM:
        {
            material_name= Material_Name::VACUUM;
            density= 0.0;
            const_eper= 1.0;
            const_mper= 1.0;
            return;
        }
    case Material_Name::CONSTANT:
        {
            material_name= Material_Name::CONSTANT;
            density= density_set;
            const_eper= const_eper_set;
            const_mper= const_mper_set;
            return;
        }
    case Material_Name::WATER:
        {
            material_name= Material_Name::WATER;
            filename= "material_water.dat";
            density= 1.0;
            break;
        }
    case Material_Name::ICE:
        {
            material_name= Material_Name::ICE;
            filename= "material_ice.dat";
            density= 0.917;
            break;
        }
    case Material_Name::BC:
        {
            material_name= Material_Name::BC;
            filename= "material_bc.dat";
            density= 1.8;
            break;
        }
    case Material_Name::MAGNETITE:
        {
            material_name= Material_Name::MAGNETITE;
            filename= "material_magnetite.dat";
            density= 5.17;
            break;
        }
    default:
            cerr << "unknown material_name" << endl;
            exit(EXIT_FAILURE);
            break;
    }
    ifstream reading_file(filename, ios::in);
    if(!reading_file){
        cerr << "file " << filename << "cannot open" << endl;
        exit(EXIT_FAILURE);
    }
    cout << "reading " << filename << "..." << endl;
    string line;
    while(getline(reading_file,line)){
        istringstream iss(line);
        string s;
        iss >> s;
        if(s == "##") continue; // comment line
        double wl0_read= stod(s);
        complex<double> eper_read;
        complex<double> mper_read;
        if(!(iss >> eper_read >> mper_read)) {
            cerr << "wrong data format in " << filename << endl;
            exit(EXIT_FAILURE);
        }
        wl0.push_back(wl0_read);
        eper.push_back(eper_read);
        mper.push_back(mper_read);
    }
}

complex<double> Material::get_eper(double free_space_wavelength) const {
    if(material_name == Material_Name::VACUUM || material_name == Material_Name::CONSTANT) return const_eper;
    if(wl0.size()<1) {
        cerr << "cannot get eper: material data is empty " << endl;
        exit(EXIT_FAILURE);
    }
    if(free_space_wavelength < wl0[0] || free_space_wavelength > wl0[wl0.size()-1]){
        cerr << "eper is unknown at this wavelength" << endl;
        exit(EXIT_FAILURE);
    }
    complex<double> eper_interp;
    for(int i= 0; i< wl0.size()-1; ++i){
        if(wl0[i] <= free_space_wavelength && free_space_wavelength <= wl0[i+1]){
            eper_interp= eper[i]+(eper[i+1]-eper[i])/(wl0[i+1]-wl0[i])*(free_space_wavelength-wl0[i]);
            return eper_interp;
        }
    }
}

complex<double> Material::get_mper(double free_space_wavelength) const {
    if(material_name == Material_Name::VACUUM || material_name == Material_Name::CONSTANT) return const_mper;
    if(wl0.size()<1) {
        cerr << "cannot get mper: material data is empty " << endl;
        exit(EXIT_FAILURE);
    }
    if(free_space_wavelength < wl0[0] || free_space_wavelength > wl0[wl0.size()-1]){
        cerr << "mper is unknown at this wavelength" << endl;
        exit(EXIT_FAILURE);
    }
    complex<double> mper_interp;
    for(int i= 0; i< wl0.size()-1; ++i){
        if(wl0[i] <= free_space_wavelength && free_space_wavelength <= wl0[i+1]){
            mper_interp= mper[i]+(mper[i+1]-mper[i])/(wl0[i+1]-wl0[i])*(free_space_wavelength-wl0[i]);
            return mper_interp;
        }
    }
}

void Material::show_info() const {
    cout << "-------- Material Info ---------" << endl;
    //cout << "Material_Name == " << material_name << endl;
    if(wl0.size()>0){
        cout << "wl0  eper  mper" << endl;
        for(int i= 0; i < wl0.size(); ++i){
            cout << wl0[i] << "  " << eper[i] << "  " << mper[i] << endl;
        }
    }else{
        cout << "const_eper const_mper " << endl;
        cout << const_eper << "  " << const_mper << endl;
    }
    cout << "-------------------------------" << endl << endl;
}
