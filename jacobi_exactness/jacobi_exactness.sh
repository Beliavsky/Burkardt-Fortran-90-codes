#! /bin/bash
#
gfortran -c -Wall jacobi_exactness.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran jacobi_exactness.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm jacobi_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/jacobi_exactness
#
echo "Normal end of execution."
