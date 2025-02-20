#! /bin/bash
#
$HOME/bin/triangulation_orient p15 > triangulation_orient_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
