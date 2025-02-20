#! /bin/bash
#
~/bin/reformat ragged.txt smooth4.txt 4
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
~/bin/reformat ragged.txt smooth7.txt 7
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
