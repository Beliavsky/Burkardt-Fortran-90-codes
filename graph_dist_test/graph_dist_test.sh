#! /bin/bash
#
gfortran -c -Wall graph_dist_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran graph_dist_test.o $HOME/lib/graph_dist.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm graph_dist_test.o
#
mv a.out graph_dist_test
./graph_dist_test > graph_dist_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm graph_dist_test
#
echo "Normal end of execution."
