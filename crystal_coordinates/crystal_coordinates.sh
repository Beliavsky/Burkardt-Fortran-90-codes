#! /bin/bash
#
gfortran -c -Wall crystal_coordinates.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran crystal_coordinates.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm crystal_coordinates.o
#
chmod ugo+x a.out
mv a.out ~/bin/crystal_coordinates
#
echo "Normal end of execution."
