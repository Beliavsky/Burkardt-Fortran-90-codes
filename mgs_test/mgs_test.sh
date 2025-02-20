#! /bin/bash
#
gfortran -c -Wall mgs_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o mgs_test mgs_test.o $HOME/lib/mgs.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm mgs_test.o
#
./mgs_test > mgs_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm mgs_test
#
echo "Normal end of execution."
