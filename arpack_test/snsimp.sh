#! /bin/bash
#
gfortran -c -Wall -I $HOME/include snsimp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran snsimp.o $HOME/lib/arpack.o -llapack -lblas
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm snsimp.o
#
mv a.out snsimp
./snsimp > snsimp.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm snsimp
#
echo "Normal end of execution."
