#! /bin/bash
#
gfortran -c -Wall prism_witherden_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o prism_witherden_rule_test prism_witherden_rule_test.o $HOME/lib/prism_witherden_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm prism_witherden_rule_test.o
#
./prism_witherden_rule_test > prism_witherden_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm prism_witherden_rule_test
#
echo "Normal end of execution."
