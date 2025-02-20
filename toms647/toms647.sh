#! /bin/bash
#
gfortran -c -Wall toms647.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms647.o ~/lib/toms647.o
#
echo "Normal end of execution."
