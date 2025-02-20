#! /bin/bash
#
gfortran -c -Wall ss_gd_align_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o ss_gd_align_test ss_gd_align_test.o $HOME/lib/ss_gd_align.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ss_gd_align_test.o
#
./ss_gd_align_test > ss_gd_align_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ss_gd_align_test
#
echo "Normal end of execution."
