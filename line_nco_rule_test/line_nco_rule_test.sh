#! /bin/bash
#
gfortran -c -Wall line_nco_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran line_nco_rule_test.o $HOME/lib/line_nco_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm line_nco_rule_test.o
#
mv a.out line_nco_rule_test
./line_nco_rule_test > line_nco_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm line_nco_rule_test
#
echo "Normal end of execution."
