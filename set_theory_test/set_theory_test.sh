#! /bin/bash
#
gfortran -c -Wall set_theory_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o set_theory_test set_theory_test.o $HOME/lib/set_theory.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm set_theory_test.o
#
./set_theory_test > set_theory_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm set_theory_test
#
echo "Normal end of execution."
