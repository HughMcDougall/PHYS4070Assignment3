//
// Created by hughm on 16/05/2023.
//
//--------------------
#include <cassert>
#include <iostream>
#include <fstream>

#include "matrix.hpp"
#include "matrix_complex.hpp"
#include "complex_vector_utils.hpp"
#include "vector_utils.hpp"
#include "dot_and_convert.hpp"
#include "physical_properties.hpp"
#include "LP_solvers.hpp"

//--------------------
using namespace matrix;
using namespace matrix_complex;

int main(int argc, char *argv[]) {
    /// Calculates time evolution of the ground state of at one 'g' value after quenching to another

    /// INPUTS
    /// Name    Type    Default             Desc
    /// tmax    double  5.0                 Maximum time of simultion
    /// dt      double  0.01                Time step size in simulation
    /// g_1     double  0.0                 Pre-quench g value to prepare in ground state
    /// g_2     double  4.0                 Post-quench g value for evolving ground state of g_1
    /// N       int     8                   Number of atoms in chain
    /// outurl  str     ./results/quenchsim location / name to save outputs to
    //-------------------------
    // Default params
    int N = 8;
    double g_1 = 0.0;
    double g_2 = 4.0;

    double tmax = 5.0;
    double dt = 0.01;
    std::string outurl = "./results/quenchsim";

    // Command line Inputs
    {
        if (argc > 1) { tmax  = std::stod(argv[1]); }
        if (argc > 2) { dt    = std::stod(argv[2]); }
        if (argc > 3) { g_1   = std::stoi(argv[3]); }
        if (argc > 4) { g_2   = std::stod(argv[4]); }
        if (argc > 5) { N     = std::stoi(argv[5]); }
        if (argc > 6) { outurl= argv[6];}
    }
    //-------------------------
    // Interpet & check inputs

    assert(tmax>0 && "max time must be > 0");
    assert(dt>0 && dt<tmax/2 && "Time step must  0 < dt <tmax");
    assert(N>1 && "No. atoms must be 2 or greater");

    int n_t = tmax / dt + 1;
    rvec T = make_rgrid(0,tmax,n_t);
    //-------------------------
    // Get ground state for pre-quench

    std::cout << "Getting ground state for g = " << g_1 << "\n";

    sqmatrix H = hamiltonian(N,g_1);

    lapack::MatrixAndVector solutions = lapack::LP_Eig_A(H);

    rvec ground_state = rvec(pow(2,N));
    for (int i=0; i<ground_state.size(); i++){
        ground_state[i] = solutions.matrix.at(0,i);
    }

    //-------------------------
    // Calculate new energy states
    std::cout << "Getting hamiltonian for g = " << g_2 << "\n";
    H = hamiltonian(N,g_2);
    solutions = lapack::LP_Eig_A(H);

    //-------------------------
    // Get U and D matrices
    std::cout << "Setting up time series...\n";

    // UDU decomposition
    sqmatrix_complex U = real_to_complex(solutions.matrix).transpose(); // Need to flip so that eigenvectors are columns, not rows
    sqmatrix_complex UT = U.transpose();
    sqmatrix_complex D(pow(2,N)); // Diagonal matrix D_ii = e^(-itE)

    // Time evolution variables
    double t = 0;
    cvec psi0 = real_to_complex(ground_state);
    cvec psi;
    sqmatrix_complex time_evolution = matrix_complex::eye(ground_state.size());

    // Observables
    sqmatrix_complex SZ = real_to_complex(S_z(N));
    sqmatrix_complex SX = real_to_complex(S_x(N));
    sqmatrix_complex CXX= real_to_complex(C_xx(N));
    double sz, sx, cxx;

    //-------------------------
    // Perform time integration

    std::ofstream outfile_real(outurl+"-real.dat");
    std::ofstream outfile_imag(outurl+"-imag.dat");

    std::ofstream outfile_Sz(outurl+"-Sz.dat");
    std::ofstream outfile_Sx(outurl+"-Sx.dat");
    std::ofstream outfile_Cxx(outurl+"-Cxx.dat");


    printv(real_d(psi0),"\t",outfile_real);
    printv(imag_d(psi0),"\t",outfile_imag);

    std::cout << "Starting time series \n";
    for(int i=0;i<n_t;i++){
        if (i%10==0){std::cout<<"\t Step "<<i<<" of "<<n_t<<"\n";}

        t = T[i];

        //Time-Evolved State
        for (int i=0; i<D.N(); i++){
            D.at(i,i) = exp(solutions.vector[i] * ci *-t);
        }
        time_evolution =  U & D & UT;  //       K = U e^lam U^T
        psi = time_evolution & psi0;   //   psi(t) = K U(0)

        // Get Observables
        sz = norm(psi & (SZ  & psi));
        sx = norm(psi & (SX  & psi));
        cxx= norm(psi & (CXX & psi));

        // Save Outputs
        printv(real_d(psi),"\t",outfile_real);
        printv(imag_d(psi),"\t",outfile_imag);

        outfile_Sz << sz  <<"\n";
        outfile_Sx << sx  <<"\n";
        outfile_Cxx<< cxx <<"\n";

    }
    outfile_real.close();
    outfile_imag.close();

    outfile_Sz.close();
    outfile_Sx.close();
    outfile_Cxx.close();

    // Save Time axis
    std::ofstream outfile_T(outurl+"-T.dat");
    printv(T, "\t", outfile_T); // Save time
    outfile_T.close();

    //-------------------------
    std::cout << "Done \n";

    return 0;
}