#! /bin/bash
#
gfortran -c -Wall glomin_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran glomin_test.o $HOME/lib/glomin.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm glomin_test.o
#
mv a.out glomin_test
./glomin_test > glomin_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm glomin_test
#
echo "Normal end of execution."
