#! /bin/bash
#
gfortran -c lapack_d.f90
#
#  Awful number of warnings!
#
#gfortran -c -Wall lapack_d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv lapack_d.o ~/lib/lapack_d.o
#
echo "Normal end of execution."
