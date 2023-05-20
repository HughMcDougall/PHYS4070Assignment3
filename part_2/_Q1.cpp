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

int main(){
    /// Presents the run-time vs atom count for the energy state calculation

    // Inputs
    std::vector<int> sizes = {2,2,3,4,5,6,7,8}; // First element as buffer due to unusual behaviour at startup
    double g = 100.0;

    // Pre-allocate time values for chrono
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Loop over different sizes
    int N;
    for(int i=0; i<sizes.size(); i++){
        N=sizes[i];
        start = std::chrono::steady_clock::now();

        //Do solution
        sqmatrix H = hamiltonian(N, g);
        lapack::MatrixAndVector solution = lapack::LP_Eig_A(H);

        end = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        if (i>0){std::cout << "Run time for N = " << N << ":\t" << duration.count() << "\n";}
    }

    return 0;
}