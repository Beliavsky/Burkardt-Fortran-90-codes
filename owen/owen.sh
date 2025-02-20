#! /bin/bash
#
gfortran -c -Wall owen.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv owen.o ~/lib/owen.o
#
echo "Normal end of execution."
