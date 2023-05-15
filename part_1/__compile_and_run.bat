:: COMPILE AND RUN .BAT FILE

echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o main main.cpp forces_and_integrators.cpp complex_vector_utils.cpp

echo "	---------------------------------"

echo "Compile complete. Running:"
echo "	---------------------------------"

:: order of inps: mode, g, out_url, u, phi, Tmax, dt, N, L, sparse
wsl ./main 2

echo "	---------------------------------"
pause