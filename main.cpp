#include <iostream>
#include <complex>
#include <vector>

#include "complex_vector_utils.hpp"

using cdouble = std::complex<double>;
using cvec = std::vector<cdouble>;

int main() {
    std::cout << "Hello, World!" << std::endl;
    cvec X(5);


    std::cout<<X[0]<<"\n";
    X+=1.0;


    return 0;
}
