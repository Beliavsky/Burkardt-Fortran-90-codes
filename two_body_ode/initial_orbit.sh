#! /bin/bash
#
gfortran -c -Wall initial_orbit.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran initial_orbit.o $HOME/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm initial_orbit.o
#
chmod ugo+x a.out
mv a.out ~/bin/initial_orbit
#
echo "Normal end of execution."
