#! /bin/bash
#
gfortran -c -Wall ps_gg_align.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ps_gg_align.o ~/lib/ps_gg_align.o
#
echo "Normal end of execution."
