#! /bin/bash
#
gfortran -c -Wall sparsekit.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparsekit.o ~/lib/sparsekit.o
#
echo "Normal end of execution."
