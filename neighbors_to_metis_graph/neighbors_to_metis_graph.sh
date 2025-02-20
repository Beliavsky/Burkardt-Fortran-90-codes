#! /bin/bash
#
gfortran -c -Wall neighbors_to_metis_graph.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran neighbors_to_metis_graph.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm neighbors_to_metis_graph.o
#
chmod ugo+x a.out
mv a.out ~/bin/neighbors_to_metis_graph
#
echo "Normal end of execution."
