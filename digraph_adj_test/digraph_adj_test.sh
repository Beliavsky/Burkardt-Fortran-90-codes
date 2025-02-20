#! /bin/bash
#
gfortran -c -g -Wall digraph_adj_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran digraph_adj_test.o $HOME/lib/digraph_adj.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm digraph_adj_test.o
#
mv a.out digraph_adj_test
./digraph_adj_test > digraph_adj_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm digraph_adj_test
#
echo "Normal end of execution."
