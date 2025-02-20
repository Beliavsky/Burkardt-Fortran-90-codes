#! /bin/bash
#
gfortran -c -Wall medit_mesh_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv medit_mesh_io.o ~/lib/medit_mesh_io.o
#
echo "Normal end of execution."
