#! /bin/bash
#
gfortran -c -Wall svd_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran svd_test.o -llapack
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm svd_test.o
#
mv a.out svd_test
mv svd_test $HOME/bin
#
echo "Normal end of execution."
