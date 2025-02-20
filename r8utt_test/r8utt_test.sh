#! /bin/bash
#
gfortran -c -Wall r8utt_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8utt_test r8utt_test.o $HOME/lib/r8utt.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8utt_test.o
#
./r8utt_test > r8utt_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8utt_test
#
echo "Normal end of execution."
