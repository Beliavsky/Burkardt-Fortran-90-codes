#! /bin/bash
#
gfortran -c -Wall matrix_chain_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv matrix_chain_brute.o ~/lib/matrix_chain_brute.o
#
echo "Normal end of execution."
