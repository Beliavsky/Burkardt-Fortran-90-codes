#! /bin/bash
#
gfortran -c -Wall quadrilateral_witherden_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o quadrilateral_witherden_rule_test quadrilateral_witherden_rule_test.o $HOME/lib/quadrilateral_witherden_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm quadrilateral_witherden_rule_test.o
#
./quadrilateral_witherden_rule_test > quadrilateral_witherden_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm quadrilateral_witherden_rule_test
#
echo "Normal end of execution."
