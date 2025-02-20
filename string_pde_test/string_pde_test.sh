#! /bin/bash
#
$HOME/bin/string_pde > string_pde_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
