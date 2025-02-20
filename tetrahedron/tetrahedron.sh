#! /bin/bash
#
gfortran -c -Wall tetrahedron.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv tetrahedron.o ~/lib/tetrahedron.o
#
echo "Normal end of execution."
