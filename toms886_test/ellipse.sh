#! /bin/bash
#
gfortran -c -Wall ellipse.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran ellipse.o $HOME/lib/toms886.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ellipse.o
#
mv a.out ellipse
./ellipse > ellipse.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ellipse
#
echo "Normal end of execution."
