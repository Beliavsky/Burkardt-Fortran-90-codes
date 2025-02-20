#! /bin/bash
#
gfortran -c -Wall graph_dist.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv graph_dist.o ~/lib/graph_dist.o
#
echo "Normal end of execution."
