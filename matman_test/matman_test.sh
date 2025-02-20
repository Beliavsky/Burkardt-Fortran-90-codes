#! /bin/bash
#
$HOME/bin/matman < input.txt > matman_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
