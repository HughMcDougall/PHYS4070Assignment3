#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>

#include "complex_vector_utils.hpp"
#include "forces_and_integrators.hpp"
#include "_defs.hpp"

cvec testfunc_loc(){
    cvec out(10);
    return out;
}

cdouble testfunc_loc_2(const cvec & X){
    cdouble out;
    out = X[0];
    return out;
}

int main(int argc, char *argv[]) {
    //-------------------------
    // Sim parameters
    double L = 20.0;    // Length & subdivisions
    int N = 128;

    double Tmax = 40.0; // Timescale & subdivisions
    double dt = 0.01;

    int sparse  = 1;    // Output Sparseness

    //-------------------------
    // Physical Properties

    // Type of wave to simulate
    int mode = 3;       // 1- cosine, 2- single soliton, 3- double soliton
    double g = 0.0;     // Attractive potential

    // Soliton Properties
    double u    = 0.0;    // Soliton Momentum
    double phi  = 0.0;    // Double soliton phase difference

    std::string out_url = "./results/sim_results";

    //------------------------------------------------------------------
    // INPUTS
    // order of inps: mode, g, out_url, u, phi, Tmax, dt, N, L, sparse
    {
        if (argc > 1) { mode    = std::stoi(argv[1]); }
        if (mode==1){g=0.0;}else{g=-1.0;}                   // Mode specific default g
        if (argc > 2) { g       = std::stod(argv[2]); }
        if (argc > 3) { out_url =               argv[3]; }
        if (argc > 4) { u       = std::stod(argv[4]); }
        if (argc > 5) { phi     = std::stod(argv[5]); }
        if (argc > 6) { Tmax    = std::stod(argv[6]); }
        if (argc > 7) { dt      = std::stod(argv[7]); }
        if (argc > 8) { N       = std::stoi(argv[8]); }
        if (argc > 9) { L       = std::stod(argv[9]); }
        if (argc >10) { sparse  = std::stoi(argv[10]); }
    }


    //========
    // Interpret and do safety checks

    double dx = L/(N-1);
    int M = Tmax / dt+1;

    assert(mode>=0 && mode <=3 && "Must select mode 1, 2 or 3");
    assert(Tmax>0.0 && "Max time must be positive");
    assert(dt>0.0 && dt<Tmax/2 && "Time step must be <Tmax/2");
    assert(L>0.0 && dt<Tmax/2 && "Width 'L' must be positive");
    assert(sparse>=1 && "output sparseness must be >=1");

    if(mode==1 && u!=0.0){std::cout<<"Warning! momentum 'u' only used for solitons (modes 2 or 3)";}
    if(mode!=3 && phi!=0.0){std::cout<<"Warning! Phase difference 3 only used for double solitons (modes 3)";}
    if(mode!=1 && g!=-1.0){std::cout<<"Warning! Soliton initial conditions are for g=-1.0 only";}

    //----------------------------
    // SET-UP

    // Time & position, for output later
    rvec X = make_grid(-L/2.0, L/2.0, N);
    rvec T = make_grid(0,Tmax,M);

    // Initial State
    cdouble x;
    cvec Y0(X.size());
    for (int k=0;k<X.size();k++){
        x = X[k];

        if (mode==1){
            // Mode 1: Cos Wave
            Y0[k]=cos(M_PI*x / L);

        } else if (mode==2){
            // Mode 1: Single Soliton
            Y0[k]=pow(2,0.5) * 1.0/cosh(x) * exp(u * ci * x);

        } else if (mode==3){
            // Mode 3: Double Soliton
            Y0[k]=pow(2,0.5) *
                    (1.0/cosh(x+L/4.0) * exp(u * ci * x)
                    + 1.0/cosh(x-L/4.0) * exp(-1.0*u * ci * x+ ci * phi)
                    );
        }
    }
    // Boundary Conditions
    Y0[0]  = 0;
    Y0[N-1]= 0;

    // Collapse differential function to state y'=f(y), for passing to RK4 integrator
    diff_func f = [dx,g](cvec Y){ return deriv_func(Y,dx,g);};

    // Save first state of wave
    std::ofstream outfile_real(out_url+"-real.dat");
    std::ofstream outfile_imag(out_url+"-imag.dat");

    //----------------------------
    // Sim Loop
    cvec Y = vcopy(Y0);
    printv(real_d(Y),"\t",outfile_real);
    printv(imag_d(Y),"\t",outfile_imag);

    for(int m=0;m<M;m++){
        //Calculate and apply euler step
        Y = runge_step(Y,f,dt);

        //Enforce boundary conditions
        Y[0]=0;
        Y[N-1]=0;

        //Write to outfile
        printv(real_d(Y),"\t",outfile_real);
        printv(imag_d(Y),"\t",outfile_imag);
    }

    //----------------------------
    outfile_real.close();
    outfile_imag.close();
    return 0;
}
