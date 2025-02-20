#! /bin/bash
#
gfortran -c -Wall kdv_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv kdv_exact.o ~/lib/kdv_exact.o
#
echo "Normal end of execution."
