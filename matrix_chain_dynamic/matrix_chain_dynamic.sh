#! /bin/bash
#
gfortran -c -Wall matrix_chain_dynamic.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv matrix_chain_dynamic.o ~/lib/matrix_chain_dynamic.o
#
echo "Normal end of execution."
