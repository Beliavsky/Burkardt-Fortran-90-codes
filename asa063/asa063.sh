#! /bin/bash
#
gfortran -c -Wall asa063.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv asa063.o ~/lib/asa063.o
#
echo "Normal end of execution."
