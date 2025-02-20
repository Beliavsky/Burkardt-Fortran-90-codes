#! /bin/bash
#
$HOME/bin/legendre_fast_rule 15 0.0 2.0 > legendre_fast_rule_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."

