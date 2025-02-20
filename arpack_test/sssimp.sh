#! /bin/bash
#
gfortran -c -Wall -I $HOME/include sssimp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sssimp.o $HOME/lib/arpack.o -llapack -lblas
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sssimp.o
#
mv a.out sssimp
./sssimp > sssimp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sssimp
#
echo "Normal end of execution."
