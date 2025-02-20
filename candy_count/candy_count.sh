#! /bin/bash
#
gfortran -c -Wall candy_count.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv candy_count.o ~/lib/candy_count.o
#
echo "Normal end of execution."
