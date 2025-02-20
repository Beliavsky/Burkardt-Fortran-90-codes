#! /bin/bash
#
gfortran -c -Wall sandia_rules_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sandia_rules_test.o $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sandia_rules_test.o
#
mv a.out sandia_rules_test
./sandia_rules_test > sandia_rules_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sandia_rules_test
#
echo "Normal end of execution."
