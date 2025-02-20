#! /bin/bash
#
gfortran -c -Wall treepack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o treepack_test treepack_test.o $HOME/lib/treepack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm treepack_test.o
#
./treepack_test > treepack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm treepack_test
#
echo "Normal end of execution."
