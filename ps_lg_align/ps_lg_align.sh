#! /bin/bash
#
gfortran -c -Wall ps_lg_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ps_lg_align.o ~/lib/ps_lg_align.o
#
echo "Normal end of execution."
