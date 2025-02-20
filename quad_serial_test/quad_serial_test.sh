#! /bin/bash
#
gfortran -c -Wall quad_serial_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o quad_serial_test quad_serial_test.o $HOME/lib/quad_serial.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm quad_serial_test.o
#
./quad_serial_test > quad_serial_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm quad_serial_test
#
echo "Normal end of execution."
