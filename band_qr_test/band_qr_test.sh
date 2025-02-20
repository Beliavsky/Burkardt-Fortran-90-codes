#! /bin/bash
#
gfortran -c -Wall band_qr_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o band_qr_test band_qr_test.o $HOME/lib/band_qr.o -llapack -lblas
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm band_qr_test.o
#
./band_qr_test > band_qr_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm band_qr_test
#
echo "Normal end of execution."
