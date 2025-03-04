#! /bin/bash
#
gfortran -c -cpp -Wall toms419_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms419_test toms419_test.o $HOME/lib/toms419.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms419_test.o
#
./toms419_test > toms419_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms419_test
#
echo "Normal end of execution."
