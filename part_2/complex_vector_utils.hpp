//
// Created by hughm on 13/03/2023.
//

#ifndef ASSIGNMENT1_COMPLEX_VECTOR_UTILS_H
#define ASSIGNMENT1_COMPLEX_VECTOR_UTILS_H

//
// Created by hughm on 13/03/2023.
//

#include <vector>
#include <functional>
#include <iostream>
#include "_defs.hpp"

//Overload vector operations to make direct products easier
//V-V Multiplication
cvec operator*=(cvec & a, const cvec & b);
cvec operator*(cvec a, const cvec & b);

//V-V Division
cvec operator/=(cvec & a, const cvec & b);
cvec operator/(cvec a, const cvec & b);

//V-V Addition
cvec operator+=(cvec & a, const cvec & b);
cvec operator+(cvec a, const cvec & b);

//V-V Subtraction
cvec operator-=(cvec & a, const cvec & b);
cvec operator-(cvec a, const cvec & b);

//V-D Multiplication
cvec operator*=(cvec & v, const double & a);
cvec operator*(cvec v, const double & a);
cvec operator*(const double & a, cvec v);

//V-D Division
cvec operator/=(cvec & v, const double & a);
cvec operator/(cvec v, const double & a);
cvec operator/(const double & a, cvec v);

//V-D Addition
cvec operator+=(cvec & v, const double & a);
cvec operator+(cvec v, const double & a);
cvec operator+(const double & a, cvec v);

//V-D Subtraction
cvec operator-=(cvec & v, const double & a);
cvec operator-(cvec v, const double & a);
cvec operator-(const double & a, cvec v);

//-----------------------------------------------
//V-C Multiplication
cvec operator*=(cvec & v, const cdouble & a);
cvec operator*(cvec v, const cdouble & a);
cvec operator*(const cdouble & a, cvec v);

//V-C Division
cvec operator/=(cvec & v, const cdouble & a);
cvec operator/(cvec v, const cdouble & a);
cvec operator/(const cdouble & a, cvec v);

//V-C Addition
cvec operator+=(cvec & v, const cdouble & a);
cvec operator+(cvec v, const cdouble & a);
cvec operator+(const cdouble & a, cvec v);

//V-C Subtraction
cvec operator-=(cvec & v, const cdouble & a);
cvec operator-(cvec v, const cdouble & a);
cvec operator-(const cdouble & a, cvec v);



//=======================================================
//Conversions & Integrations
cdouble vint(const cvec& a, double dx=0);
cdouble vsum(const cvec& a);
cvec vdiff(const cvec& a, double dx=0);

//=======================================================
// Complex specific functions
cvec conj(cvec X);
cvec real_c(cvec X);
cvec imag_c(cvec X);
rvec real_d(const cvec & X);
rvec imag_d(const cvec & X);
double vnorm2(const cvec& a);
double vnorm(const cvec& a);

//=======================================================
//Utility Functions
void printv(const cvec& a, const std::string & sep = ", ", std::ostream & targ = std::cout);
cvec vcopy(cvec a);
cvec zeros_like(cvec a);

//=======================================================
// Vector Products
cdouble operator&(const cvec & a, const cvec & b);

#endif //ASSIGNMENT1_VECTOR_UTILS_H