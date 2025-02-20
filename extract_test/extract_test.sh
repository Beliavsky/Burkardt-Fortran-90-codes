#! /bin/bash
#
~/bin/extract beta extract_test.f90 > extract_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
cat beta.f90 >> extract_test.txt
#
~/bin/extract theta extract_test.f90 >> extract_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
cat theta.f90 >> extract_test.txt
#
echo "Normal end of execution."
