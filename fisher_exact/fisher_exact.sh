#! /bin/bash
#
gfortran -c -Wall fisher_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv fisher_exact.o ~/lib/fisher_exact.o
#
echo "Normal end of execution."
