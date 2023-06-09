//
// header file for matrix class. Adapted from workshop 2 examples
//

#ifndef ASSIGNMENT1_MATRIX_H
#define ASSIGNMENT1_MATRIX_H

#pragma once
#include <cassert>
#include <vector>
#include <iostream>

//----------------------------------
using unt = std::size_t;
//----------------------------------
namespace matrix{
    using dtype = double;

    class sqmatrix{

    private:
        std::vector<dtype> _datavec;
        unt _N;
        unt _size;

    public:
        //Init
        sqmatrix(int N): _N(N), _size(N*N){_datavec.resize(_size);}

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

    };//sqmatrix
    // -----------------------------
    //Declare overload operators. Defined properly in cpp file
    sqmatrix operator+=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator+=(sqmatrix &a, const dtype &b);
    sqmatrix operator-=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator-=(sqmatrix &a, const dtype &b);
    sqmatrix operator*=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator*=(sqmatrix &a, const dtype &b);
    sqmatrix operator/=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator/=(sqmatrix &a, const dtype &b);

    sqmatrix operator+(sqmatrix a, const sqmatrix &b);
    sqmatrix operator+(sqmatrix a, const dtype &b);
    sqmatrix operator+(const dtype &b, sqmatrix a);
    sqmatrix operator-(sqmatrix a, const sqmatrix &b);
    sqmatrix operator-(sqmatrix a, const dtype &b);
    sqmatrix operator*(sqmatrix a, const sqmatrix &b);
    sqmatrix operator*(sqmatrix a, const dtype &b);
    sqmatrix operator*(const dtype &b, sqmatrix a);
    sqmatrix operator/(sqmatrix a, const sqmatrix &b);
    sqmatrix operator/(sqmatrix a, const dtype &b);

    // -----------------------------
    // Utility functions & Dot products
    sqmatrix eye(unt n);
    sqmatrix operator&(const sqmatrix & a, const sqmatrix & b);
    sqmatrix kron_prod(const sqmatrix & a, const sqmatrix & b);

} //Namespace

#endif //ASSIGNMENT1_MATRIX_H
