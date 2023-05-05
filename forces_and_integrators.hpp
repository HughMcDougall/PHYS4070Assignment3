//
// Created by hughm on 19/04/2023.
//

#ifndef ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP
#define ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP

#include <vector>
#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>
#include <complex>

#include "complex_vector_utils.hpp"

using cdouble = std::complex<double>;
using cvec = std::vector<cdouble>;

using function_1D  = std::function<cdouble(double)>;
using diff_func  = std::function<cvec(cvec)>;

cvec euler_step(cvec &X, const diff_func &f, double dt);
cvec runge_step(cvec &X, const diff_func &f, double dt);

#endif //ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP
