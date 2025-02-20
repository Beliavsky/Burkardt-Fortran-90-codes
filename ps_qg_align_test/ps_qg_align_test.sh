#! /bin/bash
#
gfortran -c -Wall ps_qg_align_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran ps_qg_align_test.o -L$HOME/lib -lps_qg_align
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ps_qg_align_test.o
#
mv a.out ps_qg_align_test
./ps_qg_align_test > ps_qg_align_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ps_qg_align_test
#
echo "Normal end of execution."
