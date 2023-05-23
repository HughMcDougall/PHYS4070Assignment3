//
// Created by hughm on 16/05/2023.
//
//--------------------
#include <cassert>
#include <iostream>
#include <fstream>
#include <chrono>

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
    /// Presents the run-time vs atom count for the energy state calculation

    /// INPUTS
    /// Name    Type    Default             Desc
    /// g       double  0.0                 Strength of spin-spin interaction
    //-------------------------
    // Default params

    // N vals to sweep over First element as buffer due to unusual behaviour at startup
    std::vector<int> sizes = {2,2,3,4,5,6,7,8};
    double g = 0.0;

    // Command line Inputs
    {
        if (argc > 1) { g    = std::stod(argv[1]); }
    }

    std::cout<< "Doing runtime benchmarks for g = " << g <<"\n";
    //-------------------------

    // Pre-allocate time values for chrono
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //-------------------------
    // Loop over different sizes
    int N;
    for(int i=0; i<sizes.size(); i++){

        // Start clock
        N = sizes[i];
        start = std::chrono::steady_clock::now();

        // Solve System
        sqmatrix H = hamiltonian(N, g);
        lapack::MatrixAndVector solution = lapack::LP_Eig_A(H);

        // Stop clock, do outputs
        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        if (i>0){std::cout << "Run time for N = " << N << ":\t" << duration.count() << "\n";}
    }

    //-------------------------

    std::cout<<"Done.\n";
    return 0;
}