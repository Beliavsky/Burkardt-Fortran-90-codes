#! /bin/bash
#
gfortran -c -Wall quadrature_weights_vandermonde_2d_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o quadrature_weights_vandermonde_2d_test quadrature_weights_vandermonde_2d_test.o $HOME/lib/quadrature_weights_vandermonde_2d.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm quadrature_weights_vandermonde_2d_test.o
#
./quadrature_weights_vandermonde_2d_test > quadrature_weights_vandermonde_2d_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm quadrature_weights_vandermonde_2d_test
#
echo "Normal end of execution."
