#! /bin/bash
#
gfortran -c -Wall gfarit.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran gfarit.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm gfarit.o
#
mv a.out ~/bin/gfarit
#
echo "Normal end of execution."
