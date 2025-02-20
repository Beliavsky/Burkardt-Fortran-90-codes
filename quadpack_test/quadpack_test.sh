#! /bin/bash
#
gfortran -c -Wall quadpack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o quadpack_test quadpack_test.o $HOME/lib/quadpack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm quadpack_test.o
#
./quadpack_test > quadpack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm quadpack_test
#
echo "Normal end of execution."
