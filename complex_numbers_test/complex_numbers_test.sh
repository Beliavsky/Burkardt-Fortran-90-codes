#! /bin/bash
#
gfortran -c -Wall complex_numbers_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran complex_numbers_test.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm complex_numbers_test.o
#
chmod ugo+x a.out
mv a.out complex_numbers_test
./complex_numbers_test > complex_numbers_test.txt
rm complex_numbers_test
#
echo "Normal end of execution."
