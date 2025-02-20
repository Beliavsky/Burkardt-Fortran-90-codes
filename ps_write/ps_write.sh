#! /bin/bash
#
gfortran -c -Wall ps_write.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ps_write.o ~/lib/ps_write.o
#
echo "Normal end of execution."
