#! /bin/bash
#
gfortran -c -Wall r8pp_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8pp_test r8pp_test.o $HOME/lib/r8pp.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8pp_test.o
#
./r8pp_test > r8pp_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8pp_test
#
echo "Normal end of execution."
