#! /bin/bash
#
gfortran -c -Wall digraph_arc_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran digraph_arc_test.o $HOME/lib/digraph_arc.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm digraph_arc_test.o
#
mv a.out digraph_arc_test
./digraph_arc_test > digraph_arc_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm digraph_arc_test
#
echo "Normal end of execution."
