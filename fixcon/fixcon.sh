#! /bin/bash
#
gfortran -c -Wall fixcon.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran fixcon.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm fixcon.o
#
chmod ugo+x a.out
mv a.out ~/bin/fixcon
#
echo "Normal end of execution."
