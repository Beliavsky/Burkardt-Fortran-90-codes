#! /bin/bash
#
gfortran -c -Wall rectangle.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran rectangle.o $HOME/lib/toms886.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rectangle.o
#
mv a.out rectangle
./rectangle > rectangle.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rectangle
#
echo "Normal end of execution."
