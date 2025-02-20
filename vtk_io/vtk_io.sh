#! /bin/bash
#
gfortran -c -Wall vtk_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv vtk_io.o ~/lib/vtk_io.o
#
echo "Normal end of execution."
