#! /bin/bash
#
gfortran -c -Wall prism_jaskowiec_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv prism_jaskowiec_rule.o ~/lib/prism_jaskowiec_rule.o
#
echo "Normal end of execution."
