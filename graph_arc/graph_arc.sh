#! /bin/bash
#
gfortran -c -Wall graph_arc.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv graph_arc.o ~/lib/graph_arc.o
#
echo "Normal end of execution."
