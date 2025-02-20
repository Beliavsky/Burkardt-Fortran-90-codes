#! /bin/bash
#
gfortran -c -Wall triangulation_mask.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_mask.f90"
  exit
fi
#
mv triangulation_mask.o ~/lib
#
echo "Partial program installed as ~/lib/triangulation_mask.o"
