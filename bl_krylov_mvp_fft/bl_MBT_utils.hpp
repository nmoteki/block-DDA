Eigen::MatrixXcd application_function(const Eigen::VectorXi& n, const int& f, Eigen::MatrixXi& n_ind, int level, const double& lf, const double& k0);
Eigen::VectorXcd BT_fft(const Eigen::VectorXi& n, const int& f, Eigen::MatrixXi n_ind, const int& level, const double& lf, const double& k0);
Eigen::MatrixXcd bl_BT_pad(const Eigen::VectorXi& n, const int& f, const Eigen::MatrixXcd& x);
Eigen::MatrixXcd bl_BT_reconstruct(const Eigen::VectorXi& n, const int& f, const Eigen::MatrixXcd& bz);
Eigen::VectorXcd MBT_fft_init(const Eigen::VectorXi& n, const int& f, const double& lf, const double& k0);
Eigen::MatrixXcd bl_MBT_fft_mvp(const Eigen::VectorXi& n,const int& f, const Eigen::VectorXcd& Au_til, const Eigen::MatrixXcd& P_hat, fftw_plan& plan_fwd, fftw_plan& plan_inv);
