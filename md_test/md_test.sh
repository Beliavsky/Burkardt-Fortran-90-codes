#! /bin/bash
#
gfortran -c -Wall md_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran md_test.o /$HOME/lib/md.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm md_test.o
#
mv a.out md_test
./md_test > md_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm md_test
#
echo "Normal end of execution."
