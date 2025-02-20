#! /bin/bash
#
gfortran -c -Wall medit_mesh_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o medit_mesh_io_test medit_mesh_io_test.o $HOME/lib/medit_mesh_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm medit_mesh_io_test.o
#
./medit_mesh_io_test > medit_mesh_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm medit_mesh_io_test
#
echo "Normal end of execution."
