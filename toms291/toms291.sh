#! /bin/bash
#
gfortran -c -Wall toms291.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms291.o ~/lib/toms291.o
#
echo "Normal end of execution."
