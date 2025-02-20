#! /bin/bash
#
gfortran -c -Wall latin_center_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o latin_center_test latin_center_test.o $HOME/lib/latin_center.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm latin_center_test.o
#
./latin_center_test > latin_center_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm latin_center_test
#
echo "Normal end of execution."
