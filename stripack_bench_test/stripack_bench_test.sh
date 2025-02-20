#! /bin/bash
#
~/bin/stripack_bench > stripack_bench_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
