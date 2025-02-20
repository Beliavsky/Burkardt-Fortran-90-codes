#! /bin/bash
#
gfortran -c -Wall latin_random_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o latin_random_test latin_random_test.o $HOME/lib/latin_random.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm latin_random_test.o
#
./latin_random_test > latin_random_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm latin_random_test
#
echo "Normal end of execution."
