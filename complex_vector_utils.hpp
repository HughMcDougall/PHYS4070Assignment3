//
// Created by hughm on 13/03/2023.
//

#ifndef ASSIGNMENT1_VECTOR_UTILS_H
#define ASSIGNMENT1_VECTOR_UTILS_H

//
// Created by hughm on 13/03/2023.
//

#include <vector>
#include <functional>

using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;



//Overload vector operations to make direct products easier
//V-V Multiplication
std::vector<double> operator*=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator*(std::vector<double> a, const std::vector<double> & b);

//V-V Division
std::vector<double> operator/=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator/(std::vector<double> a, const std::vector<double> & b);

//V-V Addition
std::vector<double> operator+=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator+(std::vector<double> a, const std::vector<double> & b);

//V-V Subtraction
std::vector<double> operator-=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator-(std::vector<double> a, const std::vector<double> & b);

//V-D Multiplication
std::vector<double> operator*=(std::vector<double> & v, const double & a);
std::vector<double> operator*(std::vector<double> v, const double & a);
std::vector<double> operator*(const double & a, std::vector<double> v);

//V-D Division
std::vector<double> operator/=(std::vector<double> & v, const double & a);
std::vector<double> operator/(std::vector<double> v, const double & a);
std::vector<double> operator/(const double & a, std::vector<double> v);

//V-D Addition
std::vector<double> operator+=(std::vector<double> & v, const double & a);
std::vector<double> operator+(std::vector<double> v, const double & a);
std::vector<double> operator+(const double & a, std::vector<double> v);

//V-D Subtraction
std::vector<double> operator-=(std::vector<double> & v, const double & a);
std::vector<double> operator-(std::vector<double> v, const double & a);
std::vector<double> operator-(const double & a, std::vector<double> v);

//=======================================================
//Conversions & Integrations
//Function to integrate over a vector
double vint(const std::vector<double>& a, double dx=0);
std::vector<double> vdiff(const std::vector<double>& a, double dx=0);
std::vector<double> make_grid(double rmin = 0.001, double rmax = 100, int n_grid = 101);
std::vector<double> vec_from_func(const function_1D& V, const std::vector<double>& rgrid);
double vsum(const std::vector<double>& a);
double vnorm(const std::vector<double>& a);

//=======================================================
//Utility Functions
void printv(const std::vector<double>& a);
std::vector<double> vcopy(std::vector<double> a);
std::vector<double> zeros_like(std::vector<double> a);


#endif //ASSIGNMENT1_VECTOR_UTILS_H
