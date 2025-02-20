#! /bin/bash
#
gfortran -c -Wall sandia_rules.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sandia_rules.o ~/lib/sandia_rules.o
#
echo "Normal end of execution."
