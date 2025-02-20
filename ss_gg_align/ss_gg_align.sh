#! /bin/bash
#
gfortran -c -Wall ss_gg_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ss_gg_align.o ~/lib/ss_gg_align.o
#
echo "Normal end of execution."
