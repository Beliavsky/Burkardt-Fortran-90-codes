#! /bin/bash
#
gfortran -c -Wall pce_ode_hermite.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv pce_ode_hermite.o ~/lib/pce_ode_hermite.o
#
echo "Normal end of execution."
