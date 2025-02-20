#! /bin/bash
#
gfortran -c -Wall quad_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran quad_rule_test.o $HOME/lib/quad_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm quad_rule_test.o
#
mv a.out quad_rule_test
./quad_rule_test > quad_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm quad_rule_test
#
echo "Normal end of execution."
