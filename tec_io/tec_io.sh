#! /bin/bash
#
gfortran -c -g -Wall tec_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv tec_io.o ~/lib/tec_io.o
#
echo "Normal end of execution."
