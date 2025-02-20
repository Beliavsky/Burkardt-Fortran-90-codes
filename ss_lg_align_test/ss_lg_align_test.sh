#! /bin/bash
#
gfortran -c -Wall ss_lg_align_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o ss_lg_align_test ss_lg_align_test.o $HOME/lib/ss_lg_align.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ss_lg_align_test.o
#
./ss_lg_align_test > ss_lg_align_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ss_lg_align_test
#
echo "Normal end of execution."
