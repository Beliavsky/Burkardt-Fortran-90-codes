#! /bin/bash
#
gfortran -c -Wall r83v_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r83v_test r83v_test.o $HOME/lib/r83v.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r83v_test.o
#
./r83v_test > r83v_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r83v_test
#
echo "Normal end of execution."
