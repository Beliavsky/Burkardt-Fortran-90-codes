#! /bin/bash
#
gfortran -c -Wall tiler_3d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tiler_3d.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm tiler_3d.o
mv a.out ~/bin/tiler_3d
#
echo "Normal end of execution."
