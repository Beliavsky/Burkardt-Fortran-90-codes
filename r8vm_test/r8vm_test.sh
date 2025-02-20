#! /bin/bash
#
gfortran -c -Wall r8vm_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8vm_test r8vm_test.o $HOME/lib/r8vm.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8vm_test.o
#
./r8vm_test > r8vm_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8vm_test
#
echo "Normal end of execution."
