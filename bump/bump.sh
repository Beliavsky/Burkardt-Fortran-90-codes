#! /bin/bash
#
gfortran -c -Wall bump.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran bump.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm bump.o
#
chmod ugo+x a.out
mv a.out ~/bin/bump
#
echo "Normal end of execution."
