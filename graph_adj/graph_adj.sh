#! /bin/bash
#
gfortran -c -Wall graph_adj.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv graph_adj.o ~/lib/graph_adj.o
#
echo "Normal end of execution."
