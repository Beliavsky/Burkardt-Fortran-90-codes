#! /bin/bash
#
gfortran -c -Wall pyramid_witherden_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv pyramid_witherden_rule.o ~/lib/pyramid_witherden_rule.o
#
echo "Normal end of execution."
