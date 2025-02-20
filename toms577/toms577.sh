#! /bin/bash
#
gfortran -c -Wall toms577.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms577.o ~/lib/toms577.o
#
echo "Normal end of execution."
