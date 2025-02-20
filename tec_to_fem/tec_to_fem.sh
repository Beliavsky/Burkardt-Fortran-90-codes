#! /bin/bash
#
gfortran -c -Wall tec_to_fem.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tec_to_fem.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm tec_to_fem.o
#
chmod ugo+x a.out
mv a.out ~/bin/tec_to_fem
#
echo "Normal end of execution."
