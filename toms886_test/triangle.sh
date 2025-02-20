#! /bin/bash
#
gfortran -c -Wall triangle.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran triangle.o $HOME/lib/toms886.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm triangle.o
#
mv a.out triangle
./triangle > triangle.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm triangle
#
echo "Normal end of execution."
