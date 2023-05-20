//
// Contains all physically meaningful functions for solving the 1D transverse ising model
//

#ifndef ASSIGNMENT_3_PHYSICAL_PROPERTIES_HPP
#define ASSIGNMENT_3_PHYSICAL_PROPERTIES_HPP

#include "vector_utils.hpp"
#include "complex_vector_utils.hpp"
#include "matrix.hpp"
#include "matrix_complex.hpp"
#include "dot_and_convert.hpp"
#include "LP_solvers.hpp"

//--------------------------------------------
// Spin Operators
sqmatrix sigma_z(int N, int m);
sqmatrix sigma_x(int N, int m);
//--------------------------------------
// Hamiltonian
sqmatrix hamiltonian(int N, double g);
//--------------------------------------
// Observables
sqmatrix S_z(int N);
sqmatrix S_x(int N);
sqmatrix C_xx(int N);

#endif //ASSIGNMENT_3_PHYSICAL_PROPERTIES_HPP
