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
        assert( a.N()==b.N() && "Tried to add two matrices of different sizes");
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
        assert( a.N()==b.N() && "Tried to subtract two matrices of different sizes");
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
        assert( a.N()==b.N() && "Tried to multiply two matrices of different sizes");
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
        assert( a.N()==b.N() && "Tried to divide two matrices of different sizes");
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

    // -----------------------------
    // Utility functions & Dot products
    sqmatrix eye(unt n){
        // Identity matrix generator
        sqmatrix out(n);
        for (unt i=0; i<n;i++){
            out.at(i,i)+=1.0;
        }
        return out;
    }

    sqmatrix operator&(const sqmatrix & a, const sqmatrix & b){
        /// Infix matrix product of two square matrices, called like AB = A & B
        assert(a.N()==b.N() && "matrices must be of same size to take matrix product");
        unt N = a.N();
        sqmatrix out(a.N());

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


    sqmatrix kron_prod(const sqmatrix & a, const sqmatrix & b){
        /// Returns the kronecker product of two real matrices
        unt N = a.N();
        unt M = b.N();

        sqmatrix out(N*M);

        for(unt i =0; i<N; i++){
            for(unt j =0; j<N; j++){

                // Skip empty cells. Will save time for sparse matrices
                if(a.at(i,j)==0){continue;}

                for(unt k =0; k<M; k++){
                    for(unt l =0; l<M; l++){
                        out.at(i*M+k, j*M+l) = a.at(i,j) * b.at(k,l);
                    }
                }

            }
        }

        return out;
    }

} // Namespace