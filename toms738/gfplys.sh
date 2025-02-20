#! /bin/bash
#
gfortran -c -Wall gfplys.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran gfplys.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm gfplys.o
#
mv a.out ~/bin/gfplys
#
echo "Normal end of execution."
