#! /bin/bash
#
gfortran -c -Wall sandia_cubature.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sandia_cubature.o ~/lib/sandia_cubature.o
#
echo "Normal end of execution."
