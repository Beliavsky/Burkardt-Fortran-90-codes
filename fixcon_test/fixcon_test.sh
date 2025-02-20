#! /bin/bash
#
~/bin/fixcon input.f > fixcon_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
