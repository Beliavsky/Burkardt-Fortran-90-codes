#! /bin/bash
#
$HOME/bin/football_rank < input.txt > football_rank_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
