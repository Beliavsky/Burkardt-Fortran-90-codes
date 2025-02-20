#! /bin/bash
#
gfortran -c -Wall quadrature_weights_vandermonde_2d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv quadrature_weights_vandermonde_2d.o ~/lib/quadrature_weights_vandermonde_2d.o
#
echo "Normal end of execution."
