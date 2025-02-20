#! /bin/bash
#
gfortran -c -Wall toms291_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms291_test toms291_test.o $HOME/lib/toms291.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms291_test.o
#
./toms291_test > toms291_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms291_test
#
echo "Normal end of execution."
