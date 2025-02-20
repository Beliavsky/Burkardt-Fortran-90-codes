#! /bin/bash
#
gfortran -c -Wall haar_transform_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran haar_transform_test.o $HOME/lib/haar_transform.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm haar_transform_test.o
#
mv a.out haar_transform_test
./haar_transform_test > haar_transform_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm haar_transform_test
#
echo "Normal end of execution."
