#! /bin/bash
#
gfortran -c -Wall root_rc.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv root_rc.o ~/lib/root_rc.o
#
echo "Normal end of execution."
