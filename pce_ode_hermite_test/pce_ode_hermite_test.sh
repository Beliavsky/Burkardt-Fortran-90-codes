#! /bin/bash
#
gfortran -c -Wall pce_ode_hermite_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o pce_ode_hermite_test pce_ode_hermite_test.o $HOME/lib/pce_ode_hermite.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pce_ode_hermite_test.o
#
./pce_ode_hermite_test > pce_ode_hermite_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm pce_ode_hermite_test
#
echo "Normal end of execution."
