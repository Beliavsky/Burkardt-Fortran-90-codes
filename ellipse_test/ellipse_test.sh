#! /bin/bash
#
gfortran -c -Wall ellipse_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran ellipse_test.o $HOME/lib/ellipse.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ellipse_test.o
#
mv a.out ellipse_test
./ellipse_test > ellipse_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ellipse_test
#
echo "Normal end of execution."
