#! /bin/bash
#
gfortran -c -Wall filon_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran filon_rule_test.o $HOME/lib/filon_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm filon_rule_test.o
#
mv a.out filon_rule_test
./filon_rule_test > filon_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm filon_rule_test
#
echo "Normal end of execution."
