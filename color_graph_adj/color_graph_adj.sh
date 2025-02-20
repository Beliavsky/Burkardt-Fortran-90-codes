#! /bin/bash
#
gfortran -c -Wall color_graph_adj.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv color_graph_adj.o ~/lib/color_graph_adj.o
#
echo "Normal end of execution."
