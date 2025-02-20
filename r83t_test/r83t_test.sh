#! /bin/bash
#
gfortran -c -Wall r83t_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r83t_test r83t_test.o $HOME/lib/r83t.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r83t_test.o
#
./r83t_test > r83t_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r83t_test
#
echo "Normal end of execution."
