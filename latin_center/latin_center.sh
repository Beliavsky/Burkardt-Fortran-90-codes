#! /bin/bash
#
gfortran -c -Wall latin_center.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv latin_center.o ~/lib/latin_center.o
#
echo "Normal end of execution."
