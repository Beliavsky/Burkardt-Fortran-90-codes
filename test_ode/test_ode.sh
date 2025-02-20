#! /bin/bash
#
gfortran -c -Wall test_ode.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv test_ode.o ~/lib/test_ode.o
#
echo "Normal end of execution."
