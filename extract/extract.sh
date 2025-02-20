#! /bin/bash
#
gfortran -c -Wall extract.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o extract extract.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm extract.o
#
mv extract ~/bin/extract
#
echo "Normal end of execution."
