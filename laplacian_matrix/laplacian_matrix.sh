#! /bin/bash
#
gfortran -c -Wall laplacian_matrix.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv laplacian_matrix.o ~/lib/laplacian_matrix.o
#
echo "Normal end of execution."
