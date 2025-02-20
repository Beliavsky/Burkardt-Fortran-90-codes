#! /bin/bash
#
if [ -f graph_code_test.txt ]; then
  rm graph_code_test.txt
fi
#
~/bin/graph_code mg_10x10.txt >  graph_code_test.txt
~/bin/graph_code mg_20x20.txt >> graph_code_test.txt
#
echo "Normal end of execution."
