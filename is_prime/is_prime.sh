#! /bin/bash
#
gfortran -c -Wall is_prime.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv is_prime.o ~/lib/is_prime.o
#
echo "Normal end of execution."
