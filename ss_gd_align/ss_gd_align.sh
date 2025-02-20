#! /bin/bash
#
gfortran -c -Wall ss_gd_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ss_gd_align.o ~/lib/ss_gd_align.o
#
echo "Normal end of execution."
