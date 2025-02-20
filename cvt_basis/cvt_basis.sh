#! /bin/bash
#
gfortran -c -Wall cvt_basis.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran cvt_basis.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm cvt_basis.o
#
chmod ugo+x a.out
mv a.out ~/bin/cvt_basis
#
echo "Normal end of execution."
