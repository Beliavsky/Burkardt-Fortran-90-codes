#!/bin/bash
#
gfortran -c -Wall hypercube_exactness.f90
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_exactness.f90"
  exit
fi
#
gfortran hypercube_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hypercube_exactness.o"
  exit
fi
rm hypercube_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/hypercube_exactness
#
echo "Executable installed as ~/bin/hypercube_exactness"
