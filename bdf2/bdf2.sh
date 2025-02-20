#! /bin/bash
#
gfortran -c -Wall bdf2.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv bdf2.o ~/lib/bdf2.o
#
echo "Normal end of execution."
