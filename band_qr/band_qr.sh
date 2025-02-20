#! /bin/bash
#
gfortran -c -Wall band_qr.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv band_qr.o ~/lib/band_qr.o
#
echo "Normal end of execution."
