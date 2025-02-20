#! /bin/bash
#
gfortran -c -Wall tsp_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tsp_brute.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm tsp_brute.o
#
mv a.out ~/bin/tsp_brute
#
echo "Normal end of execution."
