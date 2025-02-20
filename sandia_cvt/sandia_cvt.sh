#! /bin/bash
#
gfortran -c -Wall sandia_cvt.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sandia_cvt.o ~/lib/sandia_cvt.o
#
echo "Normal end of execution."
