#! /bin/bash
#
gfortran -c -Wall walsh_transform_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran walsh_transform_test.o $HOME/lib/walsh_transform.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm walsh_transform_test.o
#
mv a.out walsh_transform_test
./walsh_transform_test > walsh_transform_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm walsh_transform_test
#
echo "Normal end of execution."
