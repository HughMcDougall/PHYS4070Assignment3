//
// Class for square matrices
// HM Mar 23
//

#include "matrix_complex.hpp"
#include <iostream>
#include <cmath>

namespace matrix_complex{
//--------------------------------
//Overload operators (matrix-matrix)
//--------------------------------


//Addition
    sqmatrix_complex operator+=(sqmatrix_complex &a, const sqmatrix_complex &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) += b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix_complex operator+(sqmatrix_complex a, const sqmatrix_complex &b){ return a+=b;}

//Subtraction
    sqmatrix_complex operator-=(sqmatrix_complex &a, const sqmatrix_complex &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) -= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix_complex operator-(sqmatrix_complex a, const sqmatrix_complex &b){ return a-=b;}

//Multiplication
    sqmatrix_complex operator*=(sqmatrix_complex &a, const sqmatrix_complex &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) *= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix_complex operator*(sqmatrix_complex a, const sqmatrix_complex &b){ return a*=b;}

//Division
    sqmatrix_complex operator/=(sqmatrix_complex &a, const sqmatrix_complex &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) /= b.at(i,j);
            }
        }
        return(a);
    }
    sqmatrix_complex operator/(sqmatrix_complex a, const sqmatrix_complex &b){ return a/=b;}

//--------------------------------
//Overload operators (matrix-float)
//--------------------------------

//Addition
    sqmatrix_complex operator+=(sqmatrix_complex &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) += b;
            }
        }
        return(a);
    }
    sqmatrix_complex operator+(sqmatrix_complex a, const dtype &b){ return a+=b;}
    sqmatrix_complex operator+(const dtype &b, sqmatrix_complex a){ return a+=b;}

//Subtraction
    sqmatrix_complex operator-=(sqmatrix_complex &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) -= b;
            }
        }
        return(a);
    }
    sqmatrix_complex operator-(sqmatrix_complex a, const dtype &b){ return a-=b;}

//Multiplication
    sqmatrix_complex operator*=(sqmatrix_complex &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) *= b;
            }
        }
        return(a);
    }
    sqmatrix_complex operator*(sqmatrix_complex a, const dtype &b){ return a*=b;}
    sqmatrix_complex operator*(const dtype &b, sqmatrix_complex a){ return a*=b;}

//Division
    sqmatrix_complex operator/=(sqmatrix_complex &a, const dtype &b){
        for (unt i=0; i < a.N(); i++){
            for (unt j=0; j < a.N(); j++){
                a.at(i,j) /= b;
            }
        }
        return(a);
    }
    sqmatrix_complex operator/(sqmatrix_complex a, const dtype &b){ return a/=b;}

    // -----------------------------
    // Utility functions & Dot products
    sqmatrix_complex eye(unt n){
        // Identity matrix generator
        sqmatrix_complex out(n);

        for (unt i=0; i<n;i++){
            out.at(i,i)+=1.0;
        }

        return out;
    }

    sqmatrix_complex operator&(const sqmatrix_complex & a, const sqmatrix_complex & b){
        /// Infix matrix product of two square matrices, called like AB = A & B
        assert(a.N()==b.N() && "matrices must be of same size to take matrix product");
        unt N = a.N();
        sqmatrix_complex out(a.N());

        // Sweep over output entries
        for(unt i_out = 0; i_out<N; i_out++){
            for(unt j_out = 0; j_out<N; j_out++){

                // Calculate output entry
                for(unt k = 0; k<N; k++){
                    out.at(i_out,j_out)+=a.at(i_out,k)*b.at(k,j_out);
                }
                //---

            }
        }
        return out;
    }

    sqmatrix_complex exp(sqmatrix_complex const & X){
        // Identity matrix generator
        sqmatrix_complex out(X.N());

        for (unt i=0; i<X.N();i++){
            for (unt j=0; j<X.N();j++){
                out.at(i,j)=exp(X.at(i,j));
            }

        }

        return out;
    }

}