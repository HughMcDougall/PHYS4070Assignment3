echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o main main.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp -llapack -lblas -O3

echo "	---------------------------------"

echo "Compile complete. Running:"
echo "	---------------------------------"

wsl ./main

echo "	---------------------------------"
pause