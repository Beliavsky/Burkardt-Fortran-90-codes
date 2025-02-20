#! /bin/bash
#
gfortran -c -Wall graph_theory.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv graph_theory.o ~/lib/graph_theory.o
#
echo "Normal end of execution."
