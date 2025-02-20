#! /bin/bash
#
gfortran -c -Wall slatec.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv slatec.o ~/lib/slatec.o
#
echo "Normal end of execution."
