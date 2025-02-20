#! /bin/bash
#
gfortran -c -Wall graph_theory_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran graph_theory_test.o $HOME/lib/graph_theory.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm graph_theory_test.o
#
mv a.out graph_theory_test
./graph_theory_test > graph_theory_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm graph_theory_test
#
echo "Normal end of execution."
