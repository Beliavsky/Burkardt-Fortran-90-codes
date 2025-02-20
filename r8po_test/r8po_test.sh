#! /bin/bash
#
gfortran -c -Wall r8po_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8po_test r8po_test.o $HOME/lib/r8po.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8po_test.o
#
./r8po_test > r8po_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8po_test
#
echo "Normal end of execution."
