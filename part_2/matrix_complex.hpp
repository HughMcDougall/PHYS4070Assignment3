//
// header file for matrix class. Adapted from workshop 2 examples
//

#ifndef ASSIGNMENT1_MATRIX_COMPLEX_H
#define ASSIGNMENT1_MATRIX_COMPLEX_H

#pragma once
#include <cassert>
#include <vector>
#include <iostream>
#include <complex>

//----------------------------------
using unt = std::size_t;
//----------------------------------
namespace matrix_complex{
    using dtype = std::complex<double>;

    class sqmatrix_complex{

    private:
        std::vector<dtype> _datavec;
        unt _N;
        unt _size;

    public:
        //Init
        sqmatrix_complex(int N): _N(N), _size(N*N){_datavec.resize(_size);}

        //Gets
        dtype& at(unt i, unt j)        {return _datavec[i * _N + j];}            //For assigning
        dtype  at(unt i, unt j) const  {return _datavec[i * _N + j];}            //Const version
        unt N() const {return _N;}                                                 //Shield for _N
        unt size() const {return _size;}                                           //Shield for _size=_N*_N
        dtype *data() {return _datavec.data();}                                    //Shortcut to data array
        dtype &operator()(std::size_t i, std::size_t j) { return at(i, j); }       //Quick format of .at
        dtype operator()(std::size_t i, std::size_t j) const { return at(i, j); }  //Const Version

        //Misc Utilities
        void print(){
            for (unt i=0; i < _N; i++){
                for (unt j=0; j < _N; j++){
                    std::cout << _datavec[i*_N+j]<< "\t";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }

        //Trace
        dtype trace(){
            /// Returns the matrix's trace, the sum of its diagonal elements
            dtype out=0;
            for (int i=0;i<_N;i++){out+=_datavec[i * _N + i];}
            return(out);
        }

        //Diagonal
        std::vector<dtype> diag(){
            /// Returns the diagonal elements of the matrix in a vector
            std::vector<dtype> out(_N);
            for (int i=0;i<_N;i++){out[i]=_datavec[i * _N + i];}
            return(out);
        }

    };//sqmatrix_complex
    // -----------------------------
    //Declare overload operators. Defined properly in cpp file
    sqmatrix_complex operator+=(sqmatrix_complex &a, const sqmatrix_complex &b);
    sqmatrix_complex operator+=(sqmatrix_complex &a, const dtype &b);
    sqmatrix_complex operator-=(sqmatrix_complex &a, const sqmatrix_complex &b);
    sqmatrix_complex operator-=(sqmatrix_complex &a, const dtype &b);
    sqmatrix_complex operator*=(sqmatrix_complex &a, const sqmatrix_complex &b);
    sqmatrix_complex operator*=(sqmatrix_complex &a, const dtype &b);
    sqmatrix_complex operator/=(sqmatrix_complex &a, const sqmatrix_complex &b);
    sqmatrix_complex operator/=(sqmatrix_complex &a, const dtype &b);

    sqmatrix_complex operator+(sqmatrix_complex a, const sqmatrix_complex &b);
    sqmatrix_complex operator+(sqmatrix_complex a, const dtype &b);
    sqmatrix_complex operator+(const dtype &b, sqmatrix_complex a);
    sqmatrix_complex operator-(sqmatrix_complex a, const sqmatrix_complex &b);
    sqmatrix_complex operator-(sqmatrix_complex a, const dtype &b);
    sqmatrix_complex operator*(sqmatrix_complex a, const sqmatrix_complex &b);
    sqmatrix_complex operator*(sqmatrix_complex a, const dtype &b);
    sqmatrix_complex operator*(const dtype &b, sqmatrix_complex a);
    sqmatrix_complex operator/(sqmatrix_complex a, const sqmatrix_complex &b);
    sqmatrix_complex operator/(sqmatrix_complex a, const dtype &b);

    // -----------------------------
    // Utility functions & Dot products
    sqmatrix_complex eye(unt n);
    sqmatrix_complex exp(sqmatrix_complex const & X);
    sqmatrix_complex operator&(const sqmatrix_complex & a, const sqmatrix_complex & b);

} //Namespace

#endif //ASSIGNMENT1_MATRIX_H