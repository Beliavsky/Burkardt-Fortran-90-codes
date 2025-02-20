#!/bin/bash
#
gfortran -c triangulation_orient.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_orient.f90"
  exit
fi
#
gfortran triangulation_orient.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_orient.o"
  exit
fi
#
rm triangulation_orient.o
#
chmod ugo+x a.out
mv a.out ~/bin/triangulation_orient
#
echo "Executable installed as ~/bin/triangulation_orient"
