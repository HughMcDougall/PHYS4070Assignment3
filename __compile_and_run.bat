echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o main main.cpp forces_and_integrators.cpp complex_vector_utils.cpp

echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl main

echo "	---------------------------------"
pause