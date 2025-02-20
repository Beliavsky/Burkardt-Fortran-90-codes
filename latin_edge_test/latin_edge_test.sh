#! /bin/bash
#
gfortran -c -Wall latin_edge_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o latin_edge_test latin_edge_test.o $HOME/lib/latin_edge.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm latin_edge_test.o
#
./latin_edge_test > latin_edge_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm latin_edge_test
#
echo "Normal end of execution."
