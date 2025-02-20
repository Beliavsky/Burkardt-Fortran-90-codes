#! /bin/bash
#
gfortran -c -Wall treepack.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv treepack.o ~/lib/treepack.o
#
echo "Normal end of execution."
