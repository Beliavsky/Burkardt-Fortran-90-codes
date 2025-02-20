#! /bin/bash
#
gfortran -c -Wall pyramid_jaskowiec_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o pyramid_jaskowiec_rule_test pyramid_jaskowiec_rule_test.o $HOME/lib/pyramid_jaskowiec_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pyramid_jaskowiec_rule_test.o
#
./pyramid_jaskowiec_rule_test > pyramid_jaskowiec_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm pyramid_jaskowiec_rule_test
#
echo "Normal end of execution."
