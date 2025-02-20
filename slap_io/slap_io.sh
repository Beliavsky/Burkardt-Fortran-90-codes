#! /bin/bash
#
gfortran -c -Wall slap_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv slap_io.o ~/lib/slap_io.o
#
echo "Normal end of execution."
