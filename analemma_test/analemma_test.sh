#! /bin/bash
#
~/bin/analemma < analemma_input.txt > analemma_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
#  Call gnuplot to create the plots.
#
gnuplot analemma_commands.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
