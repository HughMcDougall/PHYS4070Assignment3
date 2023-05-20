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

int main(){
    int N = 8;
    double g = 10.0;

    sqmatrix H = hamiltonian(N,g);

    lapack::MatrixAndVector solutions = lapack::LP_Eig_A(H);

    rvec ground_state = rvec(pow(2,N));
    for (int i=0; i<ground_state.size(); i++){
        ground_state[i] = solutions.matrix.at(0,i);
    }
    rvec energies = solutions.vector;


    std::cout<<"Hamiltonian:\n";
    H.print();

    std::cout<<"Ground state:\n";
    printv(ground_state, ", ");
    std::cout<<"Energies:\n";
    printv(energies / N, ", ");

    std::cout<<"S_z:\n";
    std::cout<< (ground_state & (S_z(N) & ground_state)) <<"\n";

    std::cout<<"S_x:\n";
    std::cout<< (ground_state & (S_x(N) & ground_state)) <<"\n";

    std::cout<<"C_xx:\n";
    std::cout<< (ground_state & (C_xx(N) & ground_state)) <<"\n";

    return 0;
}