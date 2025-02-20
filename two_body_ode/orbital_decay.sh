#! /bin/bash
#
gfortran -c -Wall orbital_decay.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran orbital_decay.o $HOME/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm orbital_decay.o
#
chmod ugo+x a.out
mv a.out ~/bin/orbital_decay
#
echo "Normal end of execution."
