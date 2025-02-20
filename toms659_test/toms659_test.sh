#! /bin/bash
#
gfortran -c -Wall toms659_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms659_test toms659_test.o $HOME/lib/toms659.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms659_test.o
#
./toms659_test < input.txt > toms659_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms659_test
#
echo "Normal end of execution."
