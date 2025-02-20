#! /bin/bash
#
gfortran -c -Wall table_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran table_io_test.o $HOME/lib/table_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm table_io_test.o
#
mv a.out table_io_test
./table_io_test > table_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm table_io_test
#
echo "Normal end of execution."
