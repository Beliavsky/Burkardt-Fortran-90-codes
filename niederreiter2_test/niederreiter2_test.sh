#! /bin/bash
#
gfortran -c -Wall niederreiter2_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran niederreiter2_test.o $HOME/lib/niederreiter2.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm niederreiter2_test.o
#
mv a.out niederreiter2_test
./niederreiter2_test > niederreiter2_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm niederreiter2_test
#
echo "Normal end of execution."
