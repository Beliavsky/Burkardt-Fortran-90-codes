#! /bin/bash
#
gfortran -c -Wall football_scores_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o football_scores_test football_scores_test.o $HOME/lib/football_scores.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm football_scores_test.o
#
./football_scores_test > football_scores_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm football_scores_test
#
echo "Normal end of execution."
