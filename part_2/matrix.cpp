//
// Class for square matrices
// HM Mar 23
//

#include "matrix.hpp"
#include <iostream>

namespace matrix{
//--------------------------------
//Overload operators (matrix-matrix)
//--------------------------------


//Addition
    sqmatrix operator+=(sqmatrix &a, const sqmatrix &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) += b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix operator+(sqmatrix a, const sqmatrix &b){ return a+=b;}

//Subtraction
    sqmatrix operator-=(sqmatrix &a, const sqmatrix &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) -= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix operator-(sqmatrix a, const sqmatrix &b){ return a-=b;}

//Multiplication
    sqmatrix operator*=(sqmatrix &a, const sqmatrix &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) *= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix operator*(sqmatrix a, const sqmatrix &b){ return a*=b;}

//Division
    sqmatrix operator/=(sqmatrix &a, const sqmatrix &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) /= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix operator/(sqmatrix a, const sqmatrix &b){ return a/=b;}

//--------------------------------
//Overload operators (matrix-float)
//--------------------------------

//Addition
    sqmatrix operator+=(sqmatrix &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) += b;
            }
        }
        return(a);
    }
    sqmatrix operator+(sqmatrix a, const dtype &b){ return a+=b;}
    sqmatrix operator+(const dtype &b, sqmatrix a){ return a+=b;}

//Subtraction
    sqmatrix operator-=(sqmatrix &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) -= b;
            }
        }
        return(a);
    }
    sqmatrix operator-(sqmatrix a, const dtype &b){ return a-=b;}

//Multiplication
    sqmatrix operator*=(sqmatrix &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) *= b;
            }
        }
        return(a);
    }
    sqmatrix operator*(sqmatrix a, const dtype &b){ return a*=b;}
    sqmatrix operator*(const dtype &b, sqmatrix a){ return a*=b;}

//Division
    sqmatrix operator/=(sqmatrix &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) /= b;
            }
        }
        return(a);
    }
    sqmatrix operator/(sqmatrix a, const dtype &b){ return a/=b;}

}