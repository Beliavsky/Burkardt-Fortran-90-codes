#! /bin/bash
#
gfortran -c -Wall test_con.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv test_con.o ~/lib/test_con.o
#
echo "Normal end of execution."
