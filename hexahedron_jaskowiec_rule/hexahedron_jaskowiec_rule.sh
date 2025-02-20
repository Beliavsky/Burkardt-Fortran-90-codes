#! /bin/bash
#
gfortran -c -Wall hexahedron_jaskowiec_rule.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hexahedron_jaskowiec_rule.o ~/lib/hexahedron_jaskowiec_rule.o
#
echo "Normal end of execution."
