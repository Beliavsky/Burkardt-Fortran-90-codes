#! /bin/bash
#
gfortran -c -Wall wordsnake.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv wordsnake.o ~/lib/wordsnake.o
#
echo "Normal end of execution."
