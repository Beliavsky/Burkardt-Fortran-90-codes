#! /bin/bash
#
gfortran -c -Wall graph_arc_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran graph_arc_test.o $HOME/lib/graph_arc.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm graph_arc_test.o
#
mv a.out graph_arc_test
./graph_arc_test > graph_arc_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm graph_arc_test
#
echo "Normal end of execution."
