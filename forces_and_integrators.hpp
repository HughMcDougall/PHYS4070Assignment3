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

#include "_defs.hpp"
#include "complex_vector_utils.hpp"

cvec deriv_func(const cvec & X, const double & dx, const double & g = 0.0);

cvec euler_step(cvec &X, const diff_func &f, double dt);
cvec runge_step(cvec &X, const diff_func &f, double dt);

#endif //ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP
