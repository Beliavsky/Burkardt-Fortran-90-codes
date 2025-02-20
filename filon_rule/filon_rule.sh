#! /bin/bash
#
gfortran -c -Wall filon_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv filon_rule.o ~/lib/filon_rule.o
#
echo "Normal end of execution."
