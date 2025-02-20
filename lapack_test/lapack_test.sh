#! /bin/bash
#
gfortran -c -Wall lapack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran lapack_test.o -llapack
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm lapack_test.o
#
mv a.out lapack_test
./lapack_test > lapack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm lapack_test
#
echo "Normal end of execution."
