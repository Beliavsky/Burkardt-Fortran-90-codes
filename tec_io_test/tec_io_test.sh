#! /bin/bash
#
gfortran -c -Wall tec_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tec_io_test.o -L$HOME/lib -ltec_io
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm tec_io_test.o
#
mv a.out tec_io_test
./tec_io_test > tec_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm tec_io_test
#
echo "Normal end of execution."
