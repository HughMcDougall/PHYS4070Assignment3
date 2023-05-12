//
// Grab-bag of general utility functions for working with std::vectors
// Over-loads basic arithmetic (+, -, *, /) for vectors and doubles / ints
// HM Mar 23
// Update: added 'vsum' and 'vnorm' functions for assignment 2 pt1, 19/4
//

#include <iostream>
#include <vector>
#include <cassert>
#include <functional>
#include <cmath>
#include <complex>

//=======================================================
#include "_defs.hpp"

//=======================================================
//Overload vector operations to make direct products easier
//V-V Multiplication
cvec operator*=(cvec & a, const cvec & b){
    assert(a.size()==b.size() && "Tried to multiply two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]*=b[i];
    }
    return(a);
}
cvec operator*(cvec a, const cvec & b){return a*=b;}

//V-V Division
cvec operator/=(cvec & a, const cvec & b){
    assert(a.size()==b.size() && "Tried to divide two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]/=b[i];
    }
    return(a);
}
cvec operator/(cvec a, const cvec & b){return a/=b;}

//V-V Addition
cvec operator+=(cvec & a, const cvec & b){
    assert(a.size()==b.size() && "Tried to add two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]+=b[i];
    }
    return(a);
}
cvec operator+(cvec a, const cvec & b){return a+=b;}

//V-V Subtraction
cvec operator-=(cvec & a, const cvec & b){
    assert(a.size()==b.size() && "Tried to subtract two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]-=b[i];
    }
    return(a);
}
cvec operator-(cvec a, const cvec & b){return a-=b;}

//V-D Multiplication
cvec operator*=(cvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]*=a;
    }
    return(v);
}
cvec operator*(cvec v, const double & a){return(v*=a);}
cvec operator*(const double & a, cvec v){return(v*a);}

//V-D Division
cvec operator/=(cvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]/=a;
    }
    return(v);
}
cvec operator/(cvec v, const double & a){return(v/=a);}
cvec operator/(const double & a, cvec v){return(v/a);}

//V-D Addition
cvec operator+=(cvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]+=a;
    }
    return(v);
}
cvec operator+(cvec v, const double & a){return(v+=a);}
cvec operator+(const double & a, cvec v){return(v+a);}

//V-D Subtraction
cvec operator-=(cvec & v, const double & a){return(v+=a*-1);}
cvec operator-(cvec v, const double & a){return(v-=a);}
cvec operator-(const double & a, cvec v){return(v-a);}

//------------------------------------------------------------
//V-C Multiplication
cvec operator*=(cvec & v, const cdouble & a){
    for (int i=0; i<v.size(); i++){
        v[i]*=a;
    }
    return(v);
}
cvec operator*(cvec v, const cdouble & a){return(v*=a);}
cvec operator*(const cdouble & a, cvec v){return(v*a);}

//V-C Division
cvec operator/=(cvec & v, const cdouble & a){
    for (int i=0; i<v.size(); i++){
        v[i]/=a;
    }
    return(v);
}
cvec operator/(cvec v, const cdouble & a){return(v/=a);}
cvec operator/(const cdouble & a, cvec v){return(v/a);}

//V-C Addition
cvec operator+=(cvec & v, const cdouble & a){
    for (int i=0; i<v.size(); i++){
        v[i]+=a;
    }
    return(v);
}
cvec operator+(cvec v, const cdouble & a){return(v+=a);}
cvec operator+(const cdouble & a, cvec v){return(v+a);}

//V-C Subtraction
cvec operator-=(cvec & v, const cdouble & a){return(v+=a*-1.0);}
cvec operator-(cvec v, const cdouble & a){return(v-=a);}
cvec operator-(const cdouble & a, cvec v){return(v-a);}

//=======================================================
//Conversions & Integrations
//Function to integrate over a vector
cdouble vint(const cvec& a, double dx){
    /// Integrates over a vector with constant spacing. Uses simspsons rule
    int n = a.size();
    if (n%2==0){std::cerr<<"Warning! Integrating vector with even number of points: \t" << n;}
    cdouble out = 0;

    //Peform weighted summation in keeping /w simpsons rule
    out += a.front();
    for (int i=1; i<n; i+=2){ //Odd Indices
        out+=a[i]*4.0;
    }
    for (int i=2; i<n; i+=2){ //Even Indices
        out+=a[i]*2.0;
    }
    out += a.back();

    out /= 3;

    // If differental element has been provided, use it
    if (dx!=0){out*=dx;}

    return out;
}

cvec vdiff(const cvec& a, double dx){
    /// Differentiates a vector with constant spacing

    int n = a.size();
    cvec out(n);

    out[0] = a[1]-a[0];             //forward difference
    for (int i=1; i<n-1; i++){
        out[i]=(a[i+1]-a[i-1])/2.0;       //central difference
    }
    out[n-1] = a[n-1]-a[n-2];   //backward difference

    // If differental element has been provided, use it
    if (dx!=0){out/=dx;}

    return out;
}

std::vector<double> make_grid(double rmin = 0.001, double rmax = 100, int n_grid = 101){
    /// Creates a linspace of radial grid points. Saves time on passing the same args to every function
    double dr = (rmax-rmin) / (n_grid-1);
    double r = rmin;
    std::vector<double> out(n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=r;
        r+=dr;
    }

    return out;
}

cdouble vsum(const cvec& a){
    /// Summates a complex vector
    int n = a.size();

    cdouble out = 0;
    for (int i=0; i<n; i++){ //Odd Indices
        out+=a[i];
    }
    return out;
}

cvec vec_from_func(const function_1D& V, const std::vector<double>& rgrid){
    ///Converts a 1D potential function into a std::vector for quick-swapping of potential functions
    int n_grid = rgrid.size();
    cvec out(n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=V(rgrid[i]);
    }

    return out;
}

//=======================================================
// Complex specific functions

cvec conj(cvec X){
    /// Calculate conjugate for all elements in a vector
    for (int k=0;k<X.size();k++){X[k]=conj(X[k]);}
    return X;
}

cvec real_c(cvec X){
    /// Return real component for all elements in a vector
    for (int k=0;k<X.size();k++){X[k]=real(X[k]);}
    return X;
}

cvec imag_c(cvec X){
    /// Return imaginary component for all elements in a vector
    for (int k=0;k<X.size();k++){X[k]=imag(X[k]);}
    return X;
}

rvec real_d(const cvec & X){
    /// Return real component for all elements in a vector as a vector of doubles
    rvec out(X.size());
    for (int k=0;k<X.size();k++){out[k]=(double)real(X[k]);}
    return out;
}

rvec imag_d(const cvec & X){
    /// Return imaginary component for all elements in a vector as a vector of doubles
    rvec out(X.size());
    for (int k=0;k<X.size();k++){out[k]=(double)imag(X[k]);}
    return out;
}

double vnorm2(const cvec& a){
    /// gets the square-normal of a complex vector
    double out =0;
    for (int k=0;k<a.size();k++){
        out += norm(a[k]);
    }
    return out;
}
double vnorm(const cvec& a){
    /// gets the normal of a complex vector
    return pow(vnorm2(a),0.5);
}


//=======================================================
//Utility Functions
void printv(const cvec& a, const std::string & sep, std::ostream & targ){
    /// Prints all elements in a 1D vector of complex doubles
    for (int i=0; i<a.size(); i++){
        targ<< a[i] << sep;
    }
    targ<<"\n";
}

void printv(const rvec& a, const std::string & sep, std::ostream & targ){
    /// Prints all elements in a 1D vector of doubles
    for (int i=0; i<a.size(); i++){
        targ<< a[i] << sep;
    }
    targ<<"\n";
}

cvec vcopy(cvec a){
    /// Returns a copy of vector to assign a second version of it.
    return(a);
}

cvec zeros_like(cvec a){
    /// Returns a copy of vector to assign a second version of it.
    return(a*0);
}
