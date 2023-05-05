//
// This file contains all physically meaningful functions for PHYS4070 Assignment 2 part 1
//

#include "forces_and_integrators.hpp"

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
cdouble ci(0.0,1.0);

//================================================================
// DERIVATIVE FUNCTION
cvec deriv_func(const cvec & X, const double & dx){
    // Get Derivatives

}

//================================================================
// INTEGRATORS
cvec euler_step(cvec &X, const diff_func &f, double dt){
    /// Updates a system vector using euler step
    /// For debug purposes only
    cvec dx = f(X) * dt;
    return X + dx;
}

cvec runge_step(cvec &X, const diff_func &f, double dt){
    /// Updates a system vector using the fourth order runge kutta method

    cvec k1 = f(X);
    cvec k2 = f(X + dt/2.0 * k1);
    cvec k3 = f(X + dt/2.0 * k2);
    cvec k4 = f(X + dt   * k3);

    return X + (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0) * dt;
}
