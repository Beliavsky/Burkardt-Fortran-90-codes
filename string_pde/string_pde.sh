#! /bin/bash
#
gfortran -c -Wall string_pde.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran string_pde.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm string_pde.o
mv a.out ~/bin/string_pde
#
echo "Normal end of execution."
