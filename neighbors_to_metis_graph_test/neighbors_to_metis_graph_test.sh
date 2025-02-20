#! /bin/bash
#
$HOME/bin/neighbors_to_metis_graph cube_tetra.neighbors > neighbors_to_metis_graph_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
