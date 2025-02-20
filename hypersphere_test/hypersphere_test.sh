#! /bin/bash
#
gfortran -c -Wall hypersphere_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran hypersphere_test.o $HOME/lib/hypersphere.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hypersphere_test.o
#
mv a.out hypersphere_test
./hypersphere_test > hypersphere_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hypersphere_test
#
echo "Normal end of execution."
