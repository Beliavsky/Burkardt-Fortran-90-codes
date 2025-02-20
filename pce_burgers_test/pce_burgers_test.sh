#! /bin/bash
#
~/bin/pce_burgers > pce_burgers_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
