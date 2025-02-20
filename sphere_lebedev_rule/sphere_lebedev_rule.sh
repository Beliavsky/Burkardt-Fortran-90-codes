#! /bin/bash
#
gfortran -c -Wall sphere_lebedev_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sphere_lebedev_rule.o ~/lib/sphere_lebedev_rule.o
#
echo "Normal end of execution."
