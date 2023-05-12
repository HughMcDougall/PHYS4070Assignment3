//
// Created by hughm on 5/05/2023.
//

#ifndef ASSIGNMENT_3__DEFS_HPP
#define ASSIGNMENT_3__DEFS_HPP

#pragma once

#include <complex>
#include <vector>
#include <functional>

using cdouble = std::complex<double>;
using cvec = std::vector<cdouble>;
using rvec = std::vector<double>;

using function_1D  = std::function<cdouble(double)>;
using diff_func  = std::function<cvec(cvec)>;
inline cdouble ci(0.0,1.0);

#endif //ASSIGNMENT_3__DEFS_HPP


