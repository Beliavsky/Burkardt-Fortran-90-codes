#! /bin/bash
#
gfortran -c -Wall ss_qg_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ss_qg_align.o ~/lib/ss_qg_align.o
#
echo "Normal end of execution."
