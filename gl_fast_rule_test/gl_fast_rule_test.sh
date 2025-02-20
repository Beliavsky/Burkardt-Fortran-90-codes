#! /bin/bash
#
gfortran -c -Wall gl_fast_rule_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o gl_fast_rule_test gl_fast_rule_test.o $HOME/lib/gl_fast_rule.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm gl_fast_rule_test.o
#
./gl_fast_rule_test > gl_fast_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm gl_fast_rule_test
#
echo "Normal end of execution."
