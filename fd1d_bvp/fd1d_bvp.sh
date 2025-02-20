#! /bin/bash
#
gfortran -c -Wall fd1d_bvp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv fd1d_bvp.o ~/lib/fd1d_bvp.o
#
echo "Normal end of execution."
