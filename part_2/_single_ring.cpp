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
    /// A test file to inspect outputs for a single run

    /// INPUTS
    /// Name    Type    Default             Desc
    /// g       double  0.0                 Strength of spin-spin interaction
    //-------------------------
    // Default params

    int N = 3;
    double g = 0.0;

    // Command line Inputs
    {
        if (argc > 1) { g     = std::stod(argv[1]); }
        if (argc > 2) { N     = std::stoi(argv[2]); }
    }

    //-------------------------

    assert(N>1 && "No. atoms must be 2 or greater");

    //-------------------------
    // Solve System

    sqmatrix H = hamiltonian(N,g);

    lapack::MatrixAndVector solutions = lapack::LP_Eig_A(H);

    rvec ground_state = rvec(pow(2,N));
    for (int i=0; i<ground_state.size(); i++){
        ground_state[i] = solutions.matrix.at(0,i);
    }

    //-------------------------
    // Print Interesting outputs
    std::cout<<"Hamiltonian:\n";
    H.print();

    std::cout<<"Eigenstate Matrix:\n";
    solutions.matrix.print();

    std::cout<<"Eigenvalues / Energies:\n";
    printv(solutions.vector, "\n");

    std::cout<<"Ground State:\n";
    printv(ground_state, "\t");

    //-------------------------
    std::cout<<"Done.\n";

    return 0;
}