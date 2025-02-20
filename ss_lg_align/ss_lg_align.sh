#! /bin/bash
#
gfortran -c -Wall ss_lg_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ss_lg_align.o ~/lib/ss_lg_align.o
#
echo "Normal end of execution."
