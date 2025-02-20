#! /bin/bash
#
gfortran -c -Wall vtk_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o vtk_io_test vtk_io_test.o $HOME/lib/vtk_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm vtk_io_test.o
#
./vtk_io_test > vtk_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm vtk_io_test
#
echo "Normal end of execution."
