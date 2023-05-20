//
// Created by hughm on 16/05/2023.
//

#include "dot_and_convert.hpp"


//----------------------------------------------
using namespace matrix;
using namespace matrix_complex;

//----------------------------------------------
// Matrix Conversion
sqmatrix_complex real_to_complex(const sqmatrix & X){
    /// Conversion of real matrix to complex matrix. Changes only data type
    int N = X.N();
    sqmatrix_complex out(N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out.at(i,j) += X.at(i,j);
        }
    }

    return out;
}

sqmatrix complex_to_real(const sqmatrix_complex & X){
    /// Conversion of complex matrix to real matrix. Converts A_ij -> norm(A_ij)

    int N = X.N();
    sqmatrix out(N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out.at(i,j) = norm(X.at(i,j));
        }
    }

    return out;
}

sqmatrix_complex operator&(const sqmatrix_complex & A, const sqmatrix & B){
    return real_to_complex(B) & A;
}
sqmatrix_complex operator&(const sqmatrix & A, const sqmatrix_complex & B){return B & A;}

//----------------------------------------------
//Vector conversion / overloads
cvec operator+=(cvec & a, const rvec & b){
    assert(a.size()==b.size() && "Attempted to add rvec and cvec with different sizes");
    for(int i=0;i<a.size();i++){a[i]+=b[i];}
    return a;
}
cvec operator+(cvec a, const rvec & b){return a+=b;}
cvec operator+(const rvec & a, cvec b){return b+=a;}

cvec operator*(const rvec & a, const cdouble & x){
    cvec out(a.size());
    for (int i=0; i<a.size(); i++){out[i] = x * a[i];}
    return out;
}

cvec real_to_complex(const rvec & X){
    /// Conversion of real vector to complex vector
    int N = X.size();
    cvec out(N);
    for(int i=0; i<N; i++){
            out[i] += X[i];
    }
    return out;
}

rvec complex_to_real(const cvec & X){
    /// Conversion of complex vector to real vector. Returns as v_i = norm(x_i)
    int N = X.size();
    rvec out(N);
    for(int i=0; i<N; i++){
        out[i] = norm(X[i]);
    }
    return out;
}
//----------------------------------------------
//Matrix-Vector Products
cvec operator&(const sqmatrix_complex & M, const cvec & V){
    /// Template function for matrix multiplication of vectors, called like v = M & v

    assert(M.N()== V.size() && "Tried to take product of matrix and vector with different sizes");

    int N = M.N();
    cvec out(N);

    // Do matrix mult
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out[i] += M.at(i,j) * V[j];
        }
    }

    return out;
}

rvec operator&(const sqmatrix & M, const rvec & V){
    /// Template function for matrix multiplication of vectors, called like v = M & v

    assert(M.N()== V.size() && "Tried to take product of matrix and vector with different sizes");

    int N = M.N();
    rvec out(N);

    // Do matrix mult
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            out[i] += M.at(i,j) * V[j];
        }
    }

    return out;
}

cvec operator&(const sqmatrix_complex & M, const rvec & V){ return M & real_to_complex(V);}
cvec operator&(const sqmatrix & M, const cvec & V){ return real_to_complex(M) & V;}