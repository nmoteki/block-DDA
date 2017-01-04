std::tuple<std::vector<std::complex<double>>,std::vector<std::complex<double>>,std::vector<std::complex<double>>,std::vector<std::complex<double>>> gmie_coeff(const double& wl_0, const double& r_p, const std::complex<double>& eper_p, const std::complex<double>& mper_p, const std::complex<double>& eper_m,const std::complex<double>& mper_m);

std::tuple<double,double,double,std::vector<std::complex<double>>,std::vector<std::complex<double>>>
gmie(const double& wl_0, const double& r_p, const int& nang,
const std::complex<double>& eper_p, const std::complex<double>& mper_p,
const std::complex<double>& eper_m,const std::complex<double>& mper_m);
