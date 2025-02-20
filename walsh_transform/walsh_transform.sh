#! /bin/bash
#
gfortran -c -Wall walsh_transform.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv walsh_transform.o ~/lib/walsh_transform.o
#
echo "Normal end of execution."
