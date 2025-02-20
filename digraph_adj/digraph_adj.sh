#! /bin/bash
#
gfortran -c -g -Wall digraph_adj.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv digraph_adj.o ~/lib/digraph_adj.o
#
echo "Normal end of execution."
