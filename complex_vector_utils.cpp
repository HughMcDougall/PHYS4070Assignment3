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

using cdouble = std::complex<double>;
using cvec = std::vector<cdouble>;

using function_1D  = std::function<cdouble(double)>;
using list_of_vecs = std::vector<cvec>;

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

//Function to integrate over a vector
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

cdouble vsum(const cvec& a){
    /// Summates a vector
    int n = a.size();

    cdouble out = 0;
    for (int i=0; i<n; i++){ //Odd Indices
        out+=a[i];
    }
    return out;
}

cdouble vnorm(const cvec& a){
    /// gets the normal of a complex vector
    return pow(vsum(a * a),0.5);
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
//Utility Functions
void printv(const cvec& a){
    /// Prints all elements in a 1D vector of doubles to std::cout. Good for debugging
    for (int i=0; i<a.size(); i++){
        std::cout << a[i] << ", ";
    }
    std::cout<<"\n";
}

cvec vcopy(cvec a){
    /// Returns a copy of vector to assign a second version of it.
    return(a);
}

cvec zeros_like(cvec a){
    /// Returns a copy of vector to assign a second version of it.
    return(a*0);
}
