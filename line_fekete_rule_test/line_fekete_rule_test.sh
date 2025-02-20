#! /bin/bash
#
gfortran -c -Wall line_fekete_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran line_fekete_rule_test.o $HOME/lib/line_fekete_rule.o $HOME/lib/qr_solve.o $HOME/lib/r8lib.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm line_fekete_rule_test.o
#
mv a.out line_fekete_rule_test
./line_fekete_rule_test > line_fekete_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm line_fekete_rule_test
#
echo "Normal end of execution."
