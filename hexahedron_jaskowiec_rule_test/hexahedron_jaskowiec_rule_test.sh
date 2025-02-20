#! /bin/bash
#
gfortran -c -Wall hexahedron_jaskowiec_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o hexahedron_jaskowiec_rule_test hexahedron_jaskowiec_rule_test.o $HOME/lib/hexahedron_jaskowiec_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hexahedron_jaskowiec_rule_test.o
#
./hexahedron_jaskowiec_rule_test > hexahedron_jaskowiec_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hexahedron_jaskowiec_rule_test
#
echo "Normal end of execution."
