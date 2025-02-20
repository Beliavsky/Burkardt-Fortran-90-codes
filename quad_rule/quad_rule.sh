#! /bin/bash
#
gfortran -c -Wall quad_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv quad_rule.o ~/lib/quad_rule.o
#
echo "Normal end of execution."
