//
// Created by hughm on 16/05/2023.
//
//--------------------
#include <cassert>
#include <iostream>
#include <fstream>

#include "matrix.hpp"
#include "matrix_complex.hpp"

#include "complex_vector_utils.hpp"
#include "vector_utils.hpp"


#include "dot_and_convert.hpp"
#include "physical_properties.hpp"

#include "LP_solvers.hpp"

//--------------------
using namespace matrix;
using namespace matrix_complex;

int main(int argc, char *argv[]) {
    /// Produces ground state energies for a sweep of 'g' values for 'N' atoms

    /// INPUTS
    /// Name    Type    Default             Desc
    /// g_step  int     129                 Number of g values to sweep over
    /// g_min   double  0.0                 Min g value in sweep
    /// g_max   double  8.0                 Max g value in sweep
    /// N       int     8                   Number of atoms in chain
    /// outurl  str     ./results/energies  Location / name to save outputs to

    //-------------------------
    // Default params
    int N = 8;
    int g_step= 128 + 1;
    double g_min = 0.0;
    double g_max = 8.0;
    std::string outurl = "./results/energies";

    // Command line Inputs
    {
        if (argc > 1) { g_step= std::stoi(argv[1]); }
        if (argc > 2) { g_min = std::stod(argv[2]); }
        if (argc > 3) { g_max = std::stod(argv[3]); }
        if (argc > 4) { N     = std::stoi(argv[4]); }
        if (argc > 5) { outurl= argv[5];}
    }

    //-------------------------
    // Interpet & check inputs

    assert(g_step>0 && "number of 'g' values must be >0");
    assert(g_min!=g_max &&  "Start and end g values must be different");
    assert(N>1 && "No. atoms must be 2 or greater");

    //-------------------------
    // Set-Up

    // g values to sweep over
    rvec G = make_rgrid(g_min,g_max,g_step);
    double g;

    // Energies for storage
    rvec energies(g_step);

    //-------------------------
    // Perform Sweep
    std::cout<<"Doing energy state sweeps from g = "<< g_min << " to g = " << g_max << " with " << g_step << " steps\n";
    std::cout<<"Number of atoms = "<< N << " for grid length " << pow(2,N) << "\n";

    for (int i=0; i<g_step;i++){
        g = G[i];
        std::cout<< "\tDoing Step "<<i<<" with g = " << g <<"\n";

        // Make Hamiltonian
        sqmatrix H = hamiltonian(N,g);

        // Solve energy states
        lapack::MatrixAndVector solutions = lapack::LP_Eig_A(H);

        // Extract ground energy
        energies[i] = solutions.vector[0];
    }

    //-------------------------
    // Save outputs

    std::ofstream outfile(outurl+".dat");
    std::cout << "Calcs complete. Energies are:\n";
    std::cout << "g \t e \t e/N \n";
    for(int i=0;i<g_step;i++){
        std::cout   << G[i] << "\t" << energies[i] << "\t" << energies[i]/N << "\n";
        outfile     << G[i] << "\t" << energies[i] << "\t" << energies[i]/N << "\n";
    }
    outfile.close();

    return 0;
}