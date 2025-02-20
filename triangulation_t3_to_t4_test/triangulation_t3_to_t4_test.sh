#! /bin/bash
#
$HOME/bin/triangulation_t3_to_t4 example > triangulation_t3_to_t4_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
