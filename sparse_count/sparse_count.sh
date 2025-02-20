#! /bin/bash
#
gfortran -c -Wall sparse_count.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_count.o ~/lib/sparse_count.o
#
echo "Normal end of execution."
