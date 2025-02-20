#! /bin/bash
#
gfortran -c -Wall genin_two.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran genin_two.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm genin_two.o
#
mv a.out ~/bin/genin_two
#
echo "Normal end of execution."
