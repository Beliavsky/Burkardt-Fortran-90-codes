#! /bin/bash
#
gfortran -c -Wall pyramid_witherden_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o pyramid_witherden_rule_test pyramid_witherden_rule_test.o $HOME/lib/pyramid_witherden_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pyramid_witherden_rule_test.o
#
./pyramid_witherden_rule_test > pyramid_witherden_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm pyramid_witherden_rule_test
#
echo "Normal end of execution."
