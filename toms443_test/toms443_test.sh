#! /bin/bash
#
gfortran -c -Wall toms443_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms443_test toms443_test.o $HOME/lib/toms443.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms443_test.o
#
./toms443_test > toms443_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms443_test
#
echo "Normal end of execution."
