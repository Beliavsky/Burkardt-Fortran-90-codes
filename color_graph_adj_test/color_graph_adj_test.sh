#! /bin/bash
#
gfortran -c -Wall color_graph_adj_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran color_graph_adj_test.o $HOME/lib/color_graph_adj.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm color_graph_adj_test.o
#
mv a.out color_graph_adj_test
./color_graph_adj_test > color_graph_adj_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm color_graph_adj_test
#
echo "Normal end of execution."
