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

#include "_defs.hpp"
#include "complex_vector_utils.hpp"

//================================================================
// DERIVATIVE FUNCTION
cvec deriv_func(const cvec & X, const double & dx, const double & g){
    int N = X.size();
    double dx2 = dx*dx;
    cvec out(N);

    // Perform 2nd derivs
    for (int n=1;n<N-1;n++){
        out[n]=-1.0*(X[n-1]-2.0*X[n] + X[n+1]) / dx2;
    }

    // Pad first and last terms
    out[0]=out[1] * 0.0;
    out[N-1]=out[N-2] * 0.0;

    //Add attraction term
    out+= g * (conj(X) * X) * X;

    //Factor of i
    out*=-1.0*ci;

    return out;
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
