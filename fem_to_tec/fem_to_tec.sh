#! /bin/bash
#
gfortran -c -Wall fem_to_tec.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran fem_to_tec.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm fem_to_tec.o
#
chmod ugo+x a.out
mv a.out ~/bin/fem_to_tec
#
echo "Normal end of execution."
