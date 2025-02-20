#! /bin/bash
#
$HOME/bin/mxm 500 1000 400 > mxm_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."

