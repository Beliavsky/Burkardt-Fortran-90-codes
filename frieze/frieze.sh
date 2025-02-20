#! /bin/bash
#
gfortran -c -Wall frieze.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran frieze.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm frieze.o
#
chmod ugo+x a.out
mv a.out ~/bin/frieze
#
echo "Normal end of execution."
