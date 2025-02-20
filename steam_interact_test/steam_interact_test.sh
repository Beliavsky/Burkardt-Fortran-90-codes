#! /bin/bash
#
~/bin/steam_interact < input.txt >steam_interact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
