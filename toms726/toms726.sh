#! /bin/bash
#
gfortran -c -Wall toms726.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms726.o ~/lib/toms726.o
#
echo "Normal end of execution."
