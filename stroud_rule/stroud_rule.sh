#! /bin/bash
#
gfortran -c -Wall stroud_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv stroud_rule.o ~/lib/stroud_rule.o
#
echo "Normal end of execution."
