echo off
echo "Beginning compile"
echo "	---------------------------------"

:: wsl g++ -o _single_ring _single_ring.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3
:: wsl g++ -o _Q1 _Q1.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3
:: wsl g++ -o _Q2 _Q2.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3
wsl g++ -o _Q3 _Q3.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3

:: wsl g++ -o _testfile testfile.cpp LP_solvers.cpp vector_utils.cpp matrix.cpp complex_vector_utils.cpp matrix_complex.cpp dot_and_convert.cpp physical_properties.cpp -llapack -lblas -O3



echo "	---------------------------------"

echo "Compile complete. Running:"
echo "	---------------------------------"


:: wsl ./_Q1 
:: wsl ./_Q2 129 0 8 4 ./results/energies-4
:: wsl ./_Q2 129 0 8 5 ./results/energies-5
:: wsl ./_Q2 129 0 8 6 ./results/energies-6
:: wsl ./_Q2 129 0 8 8 ./results/energies-8
:: wsl ./_Q2 129 0 8 10 ./results/energies-10
wsl ./_Q3 
:; wsl ./_testfile



echo "	---------------------------------"
pause