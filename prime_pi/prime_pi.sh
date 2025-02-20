#! /bin/bash
#
gfortran -c -Wall prime_pi.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv prime_pi.o ~/lib/prime_pi.o
#
echo "Normal end of execution."
