#! /bin/bash
#
~/bin/tsp_brute < input.txt > tsp_brute_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
