#! /bin/bash
#
gfortran -c -Wall ps_gg_align_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o ps_gg_align_test ps_gg_align_test.o $HOME/lib/ps_gg_align.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ps_gg_align_test.o
#
./ps_gg_align_test > ps_gg_align_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ps_gg_align_test
#
echo "Normal end of execution."
