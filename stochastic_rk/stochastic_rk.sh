#! /bin/bash
#
gfortran -c -Wall stochastic_rk.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv stochastic_rk.o ~/lib/stochastic_rk.o
#
echo "Normal end of execution."
