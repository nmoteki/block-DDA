Eigen::MatrixXcd application_function(const Eigen::VectorXi& n, const int& f, Eigen::MatrixXi& n_ind, int level, const double& lf, const std::complex<double>& k);
Eigen::VectorXcd BT_fft(const Eigen::VectorXi& n, const int& f, Eigen::MatrixXi n_ind, const int& level, const double& lf, const std::complex<double>& k);
Eigen::VectorXcd BT_pad(const Eigen::VectorXi& n, const int& f, const Eigen::VectorXcd& x);
Eigen::VectorXcd BT_reconstruct(const Eigen::VectorXi& n, const int& f, const Eigen::VectorXcd& bz);
Eigen::VectorXcd MBT_fft_init(const Eigen::VectorXi& n, const int& f, const double& lf, const std::complex<double>& k,fftw_plan& plan_fwd);
Eigen::VectorXcd MBT_fft_mvp(const Eigen::VectorXi& n,const int& f, const Eigen::VectorXcd& Au_til, const Eigen::VectorXcd& p_hat, fftw_plan& plan_fwd, fftw_plan& plan_inv);
