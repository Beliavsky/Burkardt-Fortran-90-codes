#! /bin/bash
#
gfortran -c -Wall eigs_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o eigs_test eigs_test.o $HOME/lib/eigs.o -llapack
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm eigs_test.o
#
./eigs_test > eigs_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm eigs_test
#
echo "Normal end of execution."
