:: COMPILE AND RUN .BAT FILE

echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o main main.cpp forces_and_integrators.cpp complex_vector_utils.cpp

echo "	---------------------------------"

echo "Compile complete. Running:"
echo "	---------------------------------"

:: order of inps: mode, g, out_url, u, phi, Tmax, dt, N, L, sparse
wsl ./main 1 0.0 ./results/Q1-0-0
wsl ./main 1 0.1 ./results/Q1-0-1
wsl ./main 1 0.5 ./results/Q1-0-5
wsl ./main 1 2.0 ./results/Q1-2-0
wsl ./main 1 5.0 ./results/Q1-5-0

wsl ./main 2 -1.0 ./results/Q2-0 0
wsl ./main 2 -1.0 ./results/Q2--1 -1
wsl ./main 2 -1.0 ./results/Q2-1 1

wsl ./main 3 -1.0 ./results/Q3-0 0.1 0
wsl ./main 3 -1.0 ./results/Q3-0.10 0.1 0.1
wsl ./main 3 -1.0 ./results/Q3-0.25 0.1 0.25
wsl ./main 3 -1.0 ./results/Q3-0.5 0.1 0.5
wsl ./main 3 -1.0 ./results/Q3-1.0 0.1 1.0

echo "	---------------------------------"
pause