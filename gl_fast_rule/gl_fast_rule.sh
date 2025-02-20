#! /bin/bash
#
gfortran -c -Wall gl_fast_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv gl_fast_rule.o ~/lib/gl_fast_rule.o
#
echo "Normal end of execution."
