#! /bin/bash
#
gfortran -c -Wall sandia_cubature_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sandia_cubature_test sandia_cubature_test.o $HOME/lib/sandia_cubature.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sandia_cubature_test.o
#
./sandia_cubature_test > sandia_cubature_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sandia_cubature_test
#
echo "Normal end of execution."
