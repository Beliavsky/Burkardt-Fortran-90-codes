#! /bin/bash
#
gfortran -c -Wall toms660_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms660_test toms660_test.o $HOME/lib/toms660.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms660_test.o
#
./toms660_test > toms660_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms660_test
#
echo "Normal end of execution."
