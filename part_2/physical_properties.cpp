//
// Created by hughm on 20/05/2023.
//
#include "matrix.hpp"
#include "matrix_complex.hpp"

#include "vector_utils.hpp"
#include "complex_vector_utils.hpp"

#include "dot_and_convert.hpp"

#include "physical_properties.hpp"
#include "cmath"

//--------------------------------------------
using namespace matrix;
using namespace matrix_complex;

//--------------------------------------------
// Spin Operators

sqmatrix sigma_z(int N, int m){
    /// z-spin operator
    assert(m<N && "sigma_z required m in [0,N-1]");

    sqmatrix A(2);
    A.at(0,0) = 0.5;
    A.at(1,1) = -0.5;

    sqmatrix out = kron_prod(matrix::eye(pow(2,m)),
                             kron_prod(
                                     A,
                                     matrix::eye(pow(2,N-m-1))
                             )
    );

    return out;
}
sqmatrix sigma_x(int N, int m){
    /// x-spin operator

    assert(m<=N && "sigma_x required m in [0,N]");

    // Wrap-around, sigma_x(N) = sigma_x(0)
    if (m==N){ return sigma_x(N,0);}

    sqmatrix A(2);
    A.at(0,1) = 0.5;
    A.at(1,0) = 0.5;

    sqmatrix out = kron_prod(matrix::eye(pow(2,m)),
                             kron_prod(
                                     A,
                                     matrix::eye(pow(2,N-m-1))
                                       )
                             );

    return out;
}

sqmatrix _hamiltonian_int(int N){
    /// An alternate way of calculating the hamiltonian interaction term without as many matrix mults

    sqmatrix out(pow(2,N));

    sqmatrix A(2);
    A.at(0,1) = 0.5;
    A.at(1,0) = 0.5;
    sqmatrix C = kron_prod(A,A);

    out+=kron_prod(A, kron_prod(matrix::eye(pow(2,N-2)), A));
    for (int m=0; m<N-1; m++){
        out+= kron_prod(matrix::eye(pow(2,m)), kron_prod(C, matrix::eye(pow(2,N-m-2))));
    }

    return out;
}

//--------------------------------------
// Hamiltonian

sqmatrix hamiltonian(int N, double g){
    /// Calculates the hamiltonian
    sqmatrix out(pow(2,N));

    out-= _hamiltonian_int(N) * g;
    for (int m=0; m<N; m++){
        out-=sigma_z(N,m);
        //out-= g * (sigma_x(N,m) & sigma_x(N,m+1));
    }

    return out;
}

//--------------------------------------
//Observables
sqmatrix S_z(int N){
    sqmatrix out(pow(2,N));

    for (int m=0; m<N; m++){
        out+=sigma_z(N,m);
    }

    return out;
}

sqmatrix S_x(int N){
    sqmatrix out(pow(2,N));

    for (int m=0; m<N; m++){
        out+=sigma_x(N,m);
    }
    return out;
}

sqmatrix C_xx(int N){
    sqmatrix out(pow(2,N));

    for (int m=0; m<N; m++){
        for (int n=0; n<N; n++){
            if (n!=m){
                out+=sigma_x(N,m) & sigma_x(N,n) ;
            }
        }
    }
    return out;
}