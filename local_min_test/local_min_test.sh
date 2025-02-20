#! /bin/bash
#
gfortran -c -Wall local_min_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran local_min_test.o $HOME/lib/local_min.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm local_min_test.o
#
mv a.out local_min_test
./local_min_test > local_min_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm local_min_test
#
echo "Normal end of execution."
