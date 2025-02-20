#! /bin/bash
#
gfortran -c -Wall logistic_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv logistic_exact.o ~/lib/logistic_exact.o
#
echo "Normal end of execution."
