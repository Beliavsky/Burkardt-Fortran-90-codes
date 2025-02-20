#! /bin/bash
#
gfortran -c -Wall laplacian_matrix_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran laplacian_matrix_test.o $HOME/lib/laplacian_matrix.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm laplacian_matrix_test.o
#
mv a.out laplacian_matrix_test
./laplacian_matrix_test > laplacian_matrix_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm laplacian_matrix_test
#
echo "Normal end of execution."
