#! /bin/bash
#
gfortran -c -Wall porous_medium_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv porous_medium_exact.o ~/lib/porous_medium_exact.o
#
echo "Normal end of execution."
