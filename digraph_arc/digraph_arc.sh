#! /bin/bash
#
gfortran -c -Wall digraph_arc.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv digraph_arc.o ~/lib/digraph_arc.o
#
echo "Normal end of execution."
