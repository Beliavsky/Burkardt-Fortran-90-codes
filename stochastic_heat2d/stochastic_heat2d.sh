#! /bin/bash
#
gfortran -c -Wall stochastic_heat2d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv stochastic_heat2d.o ~/lib/stochastic_heat2d.o
#
echo "Normal end of execution."
