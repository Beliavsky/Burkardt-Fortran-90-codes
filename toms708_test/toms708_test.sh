#! /bin/bash
#
gfortran -c -Wall toms708_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms708_test toms708_test.o $HOME/lib/toms708.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms708_test.o
#
./toms708_test > toms708_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms708_test
#
echo "Normal end of execution."
