#! /bin/bash
#
gfortran -c -Wall mm_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv mm_io.o ~/lib/mm_io.o
#
echo "Normal end of execution."
