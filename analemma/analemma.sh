#! /bin/bash
#
gfortran -c -Wall analemma.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran analemma.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm analemma.o
#
chmod ugo+x a.out
mv a.out ~/bin/analemma
#
echo "Normal end of execution."
