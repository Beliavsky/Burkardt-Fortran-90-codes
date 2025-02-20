#! /bin/bash
#
gfortran -c -Wall md.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv md.o ~/lib/md.o
#
echo "Normal end of execution."
