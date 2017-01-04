#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include "fftw3.h"
#include <omp.h>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using namespace std::chrono;

// compile command
//g++-6 fftwtest.cpp -I/usr/local/include -L/usr/local/lib -I/Users/moteki/eigen_3_2_10 -fopenmp -lfftw3 -lfftw3_omp -lm

int main(){


    //   Typically best performance with a power of 2
 //   but could be a product of small primes
 //int  N    = 100;
 //   Compute corresponding number of complex output samples

 int  N = 134531;
 int L = 100;
 unsigned flags = FFTW_ESTIMATE;
 //unsigned flags = FFTW_MEASURE;



 cout << "omp_get_max_threads()= " << omp_get_max_threads() << endl << endl;



/****************** 1d fftw for  std::vector<compelx<double>> *******************

vector<complex<double>> a(N);
vector<complex<double>> b(N);

 for(int j=0; j<N; ++j){
     a[j]=static_cast<complex<double>>(j*0.1+1.0i);
     cout << a[j] << endl;
 }

cout << endl << endl;

fftw_plan p;
p = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&a[0]), reinterpret_cast<fftw_complex*>(&b[0]), FFTW_FORWARD, flags);



fftw_execute(p);


for(int j=0; j<N; ++j){

    cout << b[j] << endl;
}

cout << endl << endl;
cout << endl << endl;
*/
//fftw_destroy_plan(p);
//fftw_free(in);
//fftw_free(out);

/********************************************************************************



/************************ 1d fftw for  Eigen::VectorXcd *************************

VectorXcd aa(N);
VectorXcd bb(N);

fftw_plan pp;
pp = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&aa[0]), reinterpret_cast<fftw_complex*>(&bb[0]), FFTW_FORWARD, flags);
//fftw_plan ppb;
//ppb = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&bb[0]), reinterpret_cast<fftw_complex*>(&aa[0]), FFTW_BACKWARD, flags);


for(int j=0; j<N; ++j){
    aa(j)=static_cast<complex<double>>(j*0.1+1.0i);
}

cout << aa << endl;

cout << endl << endl;


fftw_execute(pp);

cout << bb << endl;

/********************************************************************************




/************************ many 1d fftws for  Eigen::MatrixXcd (each columns) using OpenMP *************************/


        {


            // N x L matrix
            MatrixXcd A(N,L);
            MatrixXcd B(N,L);


            fftw_plan p;
            p = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&(A.col(0))[0]), reinterpret_cast<fftw_complex*>(&(B.col(0))[0]), FFTW_FORWARD, flags);

            B.setZero();

            for(int j=0; j<N; ++j){
                for(int l=0; l<L; ++l){
                    A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                }
            }

            //cout << A << endl;

            omp_set_num_threads(omp_get_max_threads());

            time_point<steady_clock> t1= steady_clock::now();
            #pragma omp parallel for
            for(int l=0; l < L; ++l){

                //#pragma omp critical
                //cout << l << " " << omp_get_thread_num() << " " << endl;
                fftw_execute_dft(p,reinterpret_cast<fftw_complex*>(&(A.col(l))[0]), reinterpret_cast<fftw_complex*>(&(B.col(l))[0]));
            }
            auto t2= steady_clock::now();
            cout << "OMP fft time  (out of place) " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;

            //cout << B << endl;
            cout << endl << endl;

        }



        // in place
        {


            // N x L matrix
            MatrixXcd A(N,L);
            //MatrixXcd B(N,L);


            fftw_plan p;
            p = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&(A.col(0))[0]), reinterpret_cast<fftw_complex*>(&(A.col(0))[0]), FFTW_FORWARD, flags);

            //B.setZero();

            for(int j=0; j<N; ++j){
                for(int l=0; l<L; ++l){
                    A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                }
            }

            //cout << A << endl;

            omp_set_num_threads(omp_get_max_threads());

            time_point<steady_clock> t1= steady_clock::now();
            #pragma omp parallel for
            for(int l=0; l < L; ++l){

                //#pragma omp critical
                //cout << l << " " << omp_get_thread_num() << " " << endl;
                fftw_execute_dft(p,reinterpret_cast<fftw_complex*>(&(A.col(l))[0]), reinterpret_cast<fftw_complex*>(&(A.col(l))[0]));
            }
            auto t2= steady_clock::now();
            cout << "OMP fft time (in place) " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;

            //cout << A << endl;
            cout << endl << endl;

        }
/******************************************************************************************************************/



/************ many 1d fftws for  Eigen::Matrix<complex<double>, Dynamic,Dynamic, RowMajor> (each columns) ***************

/*
fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
fftw_complex *in, const int *inembed,
int istride, int idist,
fftw_complex *out, const int *onembed,
int ostride, int odist,
int sign, unsigned flags);
*/

/*
        {


            // Transform of each column of N x L matrix
            int rank= 1;  // 1d transform
            int n[]= {N};  // fft length
            int howmany= L;
            int idist= 1;
            int odist= 1; // distance between the first element of adjacent two fft data
            int istride= L;
            int ostride= L; // distance between two adjacent elements in each fft data
            int *inembed= n;
            int *onembed= n;

            // N x L matrix
            Matrix<complex<double>, Dynamic,Dynamic, RowMajor> A;
            Matrix<complex<double>, Dynamic,Dynamic, RowMajor> B;

            A.resize(N,L);
            B.resize(N,L);

            fftw_plan p;
            p= fftw_plan_many_dft(rank, n, howmany,
            reinterpret_cast<fftw_complex*>(&A(0,0)),
            inembed, istride, idist,
            reinterpret_cast<fftw_complex*>(&B(0,0)),
            onembed, ostride, odist,
            FFTW_FORWARD, flags);


            B.setZero();

            for(int j=0; j<N; ++j){
                for(int l=0; l<L; ++l){
                    A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                }
            }

            //cout << A << endl;
            time_point<steady_clock> t1= steady_clock::now();

            fftw_execute(p);

            auto t2= steady_clock::now();
            cout << "RowMajor many_fft time  " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;
            //cout << B << endl;
            cout << endl << endl;

        }

/******************************************************************************************************************/


/************ many 1d fftws for  Eigen::Matrix<complex<double>, Dynamic,Dynamic, ColMajor> (each columns) ***************

/*
fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
fftw_complex *in, const int *inembed,
int istride, int idist,
fftw_complex *out, const int *onembed,
int ostride, int odist,
int sign, unsigned flags);
*/

/*
        {


            // Transform of each column of N x L matrix
            int rank= 1;  // 1d transform
            int n[]= {N};  // fft length
            int howmany= L;
            int idist= N;
            int odist= N; // distance between the first element of adjacent two fft data
            int istride= 1;
            int ostride= 1; // distance between two adjacent elements in each fft data
            int *inembed= n;
            int *onembed= n;

            // N x L matrix
            /*
            Matrix<complex<double>, Dynamic,Dynamic, ColMajor> A;
            Matrix<complex<double>, Dynamic,Dynamic, ColMajor> B;

            A.resize(N,L);
            B.resize(N,L);
            */

/*
            MatrixXcd A(N,L);
            MatrixXcd B(N,L);

            fftw_plan p;
            p= fftw_plan_many_dft(rank, n, howmany,
            reinterpret_cast<fftw_complex*>(&A(0,0)),
            inembed, istride, idist,
            reinterpret_cast<fftw_complex*>(&B(0,0)),
            onembed, ostride, odist,
            FFTW_FORWARD, flags);


            B.setZero();

            for(int j=0; j<N; ++j){
                for(int l=0; l<L; ++l){
                    A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                }
            }

            //cout << A << endl;

            time_point<steady_clock> t1= steady_clock::now();

            fftw_execute(p);

            auto t2= steady_clock::now();
            cout << "ColMajor many_fft time  " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;


            //cout << B << endl;
            cout << endl << endl;
        }
        /******************************************************************************************************************/

        /************ many 1d fftws for  Eigen::Matrix<complex<double>, Dynamic,Dynamic, RowMajor> (each columns) using OpenMP ***************

        /*
        fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
        fftw_complex *in, const int *inembed,
        int istride, int idist,
        fftw_complex *out, const int *onembed,
        int ostride, int odist,
        int sign, unsigned flags);
        */

/*
                {

                    int temp= fftw_init_threads();
                    fftw_plan_with_nthreads(omp_get_max_threads());

                    // Transform of each column of N x L matrix
                    int rank= 1;  // 1d transform
                    int n[]= {N};  // fft length
                    int howmany= L;
                    int idist= 1;
                    int odist= 1; // distance between the first element of adjacent two fft data
                    int istride= L;
                    int ostride= L; // distance between two adjacent elements in each fft data
                    int *inembed= n;
                    int *onembed= n;

                    // N x L matrix
                    Matrix<complex<double>, Dynamic,Dynamic, RowMajor> A;
                    Matrix<complex<double>, Dynamic,Dynamic, RowMajor> B;

                    A.resize(N,L);
                    B.resize(N,L);

                    fftw_plan p;
                    p= fftw_plan_many_dft(rank, n, howmany,
                    reinterpret_cast<fftw_complex*>(&A(0,0)),
                    inembed, istride, idist,
                    reinterpret_cast<fftw_complex*>(&B(0,0)),
                    onembed, ostride, odist,
                    FFTW_FORWARD, flags);


                    B.setZero();

                    for(int j=0; j<N; ++j){
                        for(int l=0; l<L; ++l){
                            A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                        }
                    }

                    //cout << A << endl;
                    time_point<steady_clock> t1= steady_clock::now();

                    fftw_execute(p);

                    auto t2= steady_clock::now();
                    cout << "RowMajor many_fft _omp time  " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;
                    //cout << B << endl;
                    cout << endl << endl;

                    fftw_cleanup_threads();

                }

        /******************************************************************************************************************/


        /************ many 1d fftws for  Eigen::Matrix<complex<double>, Dynamic,Dynamic, ColMajor> (each columns)  using OpenMP ***************

        /*
        fftw_plan fftw_plan_many_dft(int rank, const int *n, int howmany,
        fftw_complex *in, const int *inembed,
        int istride, int idist,
        fftw_complex *out, const int *onembed,
        int ostride, int odist,
        int sign, unsigned flags);
        */


                {

                    int temp= fftw_init_threads();
                    fftw_plan_with_nthreads(omp_get_max_threads());

                    // Transform of each column of N x L matrix
                    int rank= 1;  // 1d transform
                    int n[]= {N};  // fft length
                    int howmany= L;
                    int idist= N;
                    int odist= N; // distance between the first element of adjacent two fft data
                    int istride= 1;
                    int ostride= 1; // distance between two adjacent elements in each fft data
                    int *inembed= n;
                    int *onembed= n;

                    // N x L matrix
                    /*
                    Matrix<complex<double>, Dynamic,Dynamic, ColMajor> A;
                    Matrix<complex<double>, Dynamic,Dynamic, ColMajor> B;

                    A.resize(N,L);
                    B.resize(N,L);
                    */


                    MatrixXcd A(N,L);
                    MatrixXcd B(N,L);

                    fftw_plan p;
                    p= fftw_plan_many_dft(rank, n, howmany,
                    reinterpret_cast<fftw_complex*>(&A(0,0)),
                    inembed, istride, idist,
                    reinterpret_cast<fftw_complex*>(&B(0,0)),
                    onembed, ostride, odist,
                    FFTW_FORWARD, flags);


                    B.setZero();

                    for(int j=0; j<N; ++j){
                        for(int l=0; l<L; ++l){
                            A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                        }
                    }

                    //cout << A << endl;

                    time_point<steady_clock> t1= steady_clock::now();

                    fftw_execute(p);

                    auto t2= steady_clock::now();
                    cout << "ColMajor many_fft _omp time (out of place) " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;


                    //cout << B << endl;
                    cout << endl << endl;

                    fftw_cleanup_threads();
                }



                // in_place
                {

                    int temp= fftw_init_threads();
                    fftw_plan_with_nthreads(omp_get_max_threads());

                    // Transform of each column of N x L matrix
                    int rank= 1;  // 1d transform
                    int n[]= {N};  // fft length
                    int howmany= L;
                    int idist= N;
                    int odist= N; // distance between the first element of adjacent two fft data
                    int istride= 1;
                    int ostride= 1; // distance between two adjacent elements in each fft data
                    int *inembed= n;
                    int *onembed= n;

                    // N x L matrix
                    /*
                    Matrix<complex<double>, Dynamic,Dynamic, ColMajor> A;
                    Matrix<complex<double>, Dynamic,Dynamic, ColMajor> B;

                    A.resize(N,L);
                    B.resize(N,L);
                    */


                    MatrixXcd A(N,L);
                    //MatrixXcd B(N,L);

                    fftw_plan p;
                    p= fftw_plan_many_dft(rank, n, howmany,
                    reinterpret_cast<fftw_complex*>(&A(0,0)),
                    inembed, istride, idist,
                    reinterpret_cast<fftw_complex*>(&A(0,0)),
                    onembed, ostride, odist,
                    FFTW_FORWARD, flags);


                    //B.setZero();

                    for(int j=0; j<N; ++j){
                        for(int l=0; l<L; ++l){
                            A(j,l)=static_cast<complex<double>>(j*(0.1+l)+l*1.0i);
                        }
                    }

                    //cout << A << endl;

                    time_point<steady_clock> t1= steady_clock::now();

                    fftw_execute(p);

                    auto t2= steady_clock::now();
                    cout << "ColMajor many_fft _omp time (in place) " << duration_cast<milliseconds>((t2 - t1)).count() << " millisecond" << endl;


                    //cout << A << endl;
                    cout << endl << endl;

                    fftw_cleanup_threads();
                }
                /******************************************************************************************************************/

}
