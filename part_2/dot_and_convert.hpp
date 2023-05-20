//
// Dot products, matrix products and real - complex conversions for moving from part 2 to part 3
//

#ifndef ASSIGNMENT_3_DOT_AND_CONVERT_HPP
#define ASSIGNMENT_3_DOT_AND_CONVERT_HPP

#include "_defs.hpp"
#include "complex_vector_utils.hpp"
#include "vector_utils.hpp"
#include "matrix.hpp"
#include "matrix_complex.hpp"

//----------------------------------------------
using namespace matrix;
using namespace matrix_complex;

//----------------------------------------------
// Matrix Conversion
sqmatrix_complex real_to_complex(const sqmatrix & X);
sqmatrix complex_to_real(const sqmatrix_complex & X);
sqmatrix_complex operator&(const sqmatrix_complex & A, const sqmatrix & B);
sqmatrix_complex operator&(const sqmatrix & A, const sqmatrix_complex & B);
//----------------------------------------------
//Vector conversion / overloads
cvec operator+=(cvec & a, const rvec & b);
cvec operator+(cvec a, const rvec & b);
cvec operator+(const rvec & a, cvec b);

cvec operator*(const rvec & a, const cdouble & x);
cvec real_to_complex(const rvec & X);
rvec complex_to_real(const cvec & X);
//----------------------------------------------
//Matrix-Vector Products
cvec operator&(const sqmatrix_complex & M, const cvec & V);
rvec operator&(const sqmatrix & M, const rvec & V);
cvec operator&(const sqmatrix_complex & M, const rvec & V);
cvec operator&(const sqmatrix & M, const cvec & V);
#endif //ASSIGNMENT_3_DOT_AND_CONVERT_HPP