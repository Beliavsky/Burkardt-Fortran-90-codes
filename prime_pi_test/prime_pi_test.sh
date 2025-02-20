#! /bin/bash
#
gfortran -c -Wall prime_pi_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o prime_pi_test prime_pi_test.o $HOME/lib/prime_pi.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm prime_pi_test.o
#
./prime_pi_test > prime_pi_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm prime_pi_test
#
echo "Normal end of execution."
