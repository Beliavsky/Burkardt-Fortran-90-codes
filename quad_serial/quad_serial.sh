#! /bin/bash
#
gfortran -c -Wall quad_serial.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv quad_serial.o ~/lib/quad_serial.o
#
echo "Normal end of execution."
