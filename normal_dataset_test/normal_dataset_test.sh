#! /bin/bash
#
~/bin/normal_dataset 2 1000 123456789 1 2 1 0 0 3 > normal_dataset_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
