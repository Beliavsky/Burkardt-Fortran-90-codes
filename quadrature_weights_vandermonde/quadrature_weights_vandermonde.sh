#! /bin/bash
#
gfortran -c -Wall quadrature_weights_vandermonde.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv quadrature_weights_vandermonde.o ~/lib/quadrature_weights_vandermonde.o
#
echo "Normal end of execution."
