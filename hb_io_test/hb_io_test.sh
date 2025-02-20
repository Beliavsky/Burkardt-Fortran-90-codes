#! /bin/bash
#
cp ~/include/hb_file_module.mod .
#
gfortran -c -Wall hb_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
rm hb_file_module.mod
#
gfortran hb_io_test.o $HOME/lib/hb_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hb_io_test.o
#
mv a.out hb_io_test
./hb_io_test > hb_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hb_io_test
#
echo "Normal end of execution."
