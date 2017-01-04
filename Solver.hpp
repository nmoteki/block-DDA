#ifndef INCLUDED_SOLVER
#define INCLUDED_SOLVER

#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include "fftw3.h"

#ifndef INCLUDED_MBT_UTILS
#include "MBT_utils.hpp"
#endif

#ifndef INCLUDED_TARGET
#include "Target.hpp"
#endif

#ifndef INCLUDED_INCIDENT_WAVE
#include "Incident_Wave.hpp"
#endif

#include "DDA_linear_solvers.hpp"


class Solver {
    Eigen::VectorXcd Au_til;
    Eigen::VectorXcd DIAG_A;
    double tol= 1e-2;
    int itermax= 300;
    int fft_length;
    fftw_plan plan_fwd;
    fftw_plan plan_inv;
    void fft_prep(const Target& target);
public:
    void show_info();
    void prepare_DDA_matrix(const Target& target, const Incident_Wave& incwave);
    void DDA_solve(Target& target);
};


#endif
