echo off
echo "Beginning compile"
echo "	---------------------------------"

:: wsl g++ -o main main.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3
wsl g++ -o _Q1 _Q1.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3

echo "	---------------------------------"

echo "Compile complete. Running:"
echo "	---------------------------------"

:: wsl ./main
wsl ./_Q1

echo "	---------------------------------"
pause