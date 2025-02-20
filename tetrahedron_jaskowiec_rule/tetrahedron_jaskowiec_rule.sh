#! /bin/bash
#
gfortran -c -Wall tetrahedron_jaskowiec_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv tetrahedron_jaskowiec_rule.o ~/lib/tetrahedron_jaskowiec_rule.o
#
echo "Normal end of execution."
