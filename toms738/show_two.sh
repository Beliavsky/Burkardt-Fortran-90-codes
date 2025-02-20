#! /bin/bash
#
gfortran -c -Wall show_two.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran show_two.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm show_two.o
#
mv a.out ~/bin/show_two
#
echo "Normal end of execution."
