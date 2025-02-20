#! /bin/bash
#
gfortran -c -Wall latin_random.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv latin_random.o ~/lib/latin_random.o
#
echo "Normal end of execution."
