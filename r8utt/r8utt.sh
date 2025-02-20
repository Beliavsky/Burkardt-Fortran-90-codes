#! /bin/bash
#
gfortran -c -Wall r8utt.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8utt.o ~/lib/r8utt.o
#
echo "Normal end of execution."
