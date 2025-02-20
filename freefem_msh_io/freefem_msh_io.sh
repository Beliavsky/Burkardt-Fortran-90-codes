#! /bin/bash
#
gfortran -c -Wall freefem_msh_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv freefem_msh_io.o ~/lib/freefem_msh_io.o
#
echo "Normal end of execution."
