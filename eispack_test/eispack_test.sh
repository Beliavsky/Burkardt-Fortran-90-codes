#! /bin/bash
#
gfortran -c -Wall eispack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran eispack_test.o $HOME/lib/eispack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm eispack_test.o
#
mv a.out eispack_test
./eispack_test > eispack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm eispack_test
#
echo "Normal end of execution."
