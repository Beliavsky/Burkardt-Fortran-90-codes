#! /bin/bash
#
gfortran -c -Wall genin.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran genin.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm genin.o
#
mv a.out ~/bin/genin
#
echo "Normal end of execution."
