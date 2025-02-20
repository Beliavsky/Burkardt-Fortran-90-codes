#! /bin/bash
#
gfortran -c -Wall tetrahedron_witherden_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv tetrahedron_witherden_rule.o ~/lib/tetrahedron_witherden_rule.o
#
echo "Normal end of execution."
