#! /bin/bash
#
gfortran -c -Wall sphere_triangle_quad.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sphere_triangle_quad.o ~/lib/sphere_triangle_quad.o
#
echo "Normal end of execution."
