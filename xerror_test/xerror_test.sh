#! /bin/bash
#
gfortran -c -Wall xerror_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o xerror_test xerror_test.o $HOME/lib/xerror.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm xerror_test.o
#
./xerror_test > xerror_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm xerror_test
#
echo "Normal end of execution."
