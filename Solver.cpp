#ifndef INCLUDED_SOLVER
#include "Solver.hpp"
#endif

using namespace Eigen;
using namespace std;

void Solver::show_info(){
    cout << "---------Solver Info---------" << endl;
    cout << "tol == " << tol << endl;
    cout << "itermax == " << itermax << endl;
    cout << "fft_length == " << fft_length << endl;
    cout << "--------------------------------------" << endl << endl;
}


void Solver::fft_prep(const Target& target){
    fft_length= target.f*target.f*(2*target.n.array()-1).prod();
    //MatrixXcd plan_buf(fft_length,1);
    VectorXcd plan_buf(fft_length);
    //plan_fwd = fftw_plan_dft_1d(fft_length, reinterpret_cast<fftw_complex*>(&(plan_buf.col(0))[0]), reinterpret_cast<fftw_complex*>(&(plan_buf.col(0))[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    //plan_inv = fftw_plan_dft_1d(fft_length, reinterpret_cast<fftw_complex*>(&(plan_buf.col(0))[0]), reinterpret_cast<fftw_complex*>(&(plan_buf.col(0))[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_fwd = fftw_plan_dft_1d(fft_length, reinterpret_cast<fftw_complex*>(&plan_buf[0]), reinterpret_cast<fftw_complex*>(&plan_buf[0]), FFTW_FORWARD, FFTW_ESTIMATE);
    plan_inv = fftw_plan_dft_1d(fft_length, reinterpret_cast<fftw_complex*>(&plan_buf[0]), reinterpret_cast<fftw_complex*>(&plan_buf[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
}


void Solver::prepare_DDA_matrix(const Target& target, const Incident_Wave& incwave){
    switch (target.f)
    {
    case 3:
        {
            DIAG_A.resize(3*target.num_element_occupy);
            for(int m= 0; m < target.num_element_occupy; ++m){
                DIAG_A.segment(3*m,3).setConstant(1.0/target.alpha_E(m));
            }
            break;
        }
    case 6:
        {
            DIAG_A.resize(6*target.num_element_occupy);
            for(int m= 0; m < target.num_element_occupy; ++m){
                DIAG_A.segment(6*m,3).setConstant(1.0/target.alpha_E(m));
                DIAG_A.segment(6*m+3,3).setConstant(1.0/target.alpha_H(m));
            }
            break;
        }
    default:
        cerr << "unknown f value" << endl;
        exit(EXIT_FAILURE);
        break;
    }

    fft_prep(target);
    Au_til = MBT_fft_init(target.n, target.f, target.lf, incwave.get_k(),plan_fwd);
}

void Solver::DDA_solve(Target& target){
    auto start_t= chrono::steady_clock::now();
    switch (target.f)
    {
    case 3:
        {
            target.PM= bl_cocg_rq_jacobi_mvp_fft(target.n, target.f, target.address, Au_til, DIAG_A, target.EHinc, tol, itermax, plan_fwd, plan_inv);
            //target.PM= bl_bicgstab_jacobi_mvp_fft(target.n, target.f, target.address, Au_til, DIAG_A, target.EHinc, tol, itermax, plan_fwd, plan_inv);
            //target.PM= bl_bicggr_jacobi_mvp_fft(target.n, target.f, target.address, Au_til, DIAG_A, target.EHinc, tol, itermax, plan_fwd, plan_inv);
            break;
        }
    case 6:
        {
            target.PM= bl_bicgstab_jacobi_mvp_fft(target.n, target.f, target.address, Au_til, DIAG_A, target.EHinc, tol, itermax, plan_fwd, plan_inv);
            //target.PM= bl_bicggr_jacobi_mvp_fft(target.n, target.f, target.address, Au_til, DIAG_A, target.EHinc, tol, itermax, plan_fwd, plan_inv);
            break;
        }
    default:
        cerr << "unknown f value" << endl;
        exit(EXIT_FAILURE);
        break;
    }
   auto end_t= chrono::steady_clock::now();
   target.DDA_solver_time= chrono::duration_cast<chrono::seconds>((end_t - start_t)).count();

}
