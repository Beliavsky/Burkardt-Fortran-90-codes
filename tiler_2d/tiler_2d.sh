#! /bin/bash
#
gfortran -c -Wall tiler_2d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tiler_2d.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm tiler_2d.o
mv a.out ~/bin/tiler_2d
#
echo "Normal end of execution."
