#! /bin/bash
#
gfortran -c -Wall r83p_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r83p_test r83p_test.o $HOME/lib/r83p.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r83p_test.o
#
./r83p_test > r83p_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r83p_test
#
echo "Normal end of execution."
