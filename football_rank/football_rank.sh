#! /bin/bash
#
gfortran -c -Wall football_rank.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran football_rank.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm football_rank.o
#
mv a.out ~/bin/football_rank
#
echo "Normal end of execution."
