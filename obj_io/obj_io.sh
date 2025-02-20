#! /bin/bash
#
gfortran -c -Wall obj_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv obj_io.o ~/lib/obj_io.o
#
echo "Normal end of execution."
