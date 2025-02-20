#! /bin/bash
#
gfortran -c -Wall grf_to_eps.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran grf_to_eps.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm grf_to_eps.o
#
chmod ugo+x a.out
mv a.out ~/bin/grf_to_eps
#
echo "Normal end of execution."
