#! /bin/bash
#
gfortran -c -Wall companion_matrix.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv companion_matrix.o ~/lib/companion_matrix.o
#
echo "Normal end of execution."
