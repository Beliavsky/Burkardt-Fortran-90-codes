#! /bin/bash
#
gfortran -c -Wall graph_adj_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran graph_adj_test.o $HOME/lib/graph_adj.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm graph_adj_test.o
#
mv a.out graph_adj_test
./graph_adj_test > graph_adj_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm graph_adj_test
#
echo "Normal end of execution."
