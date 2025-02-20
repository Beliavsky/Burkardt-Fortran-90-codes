#! /bin/bash
#
gfortran -c -Wall color_digraph_adj.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv color_digraph_adj.o ~/lib/color_digraph_adj.o
#
echo "Normal end of execution."
