#! /bin/bash
#
~/bin/tec_to_fem tiny.dat > tec_to_fem_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
