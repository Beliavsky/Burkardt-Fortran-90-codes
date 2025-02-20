#! /bin/bash
#
gfortran -c -Wall toms443.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms443.o ~/lib/toms443.o
#
echo "Normal end of execution."
