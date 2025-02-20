#! /bin/bash
#
gfortran -c -g -Wall nint_exactness_mixed.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran nint_exactness_mixed.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm nint_exactness_mixed.o
#
chmod ugo+x a.out
mv a.out ~/bin/nint_exactness_mixed
#
echo "Normal end of execution."
