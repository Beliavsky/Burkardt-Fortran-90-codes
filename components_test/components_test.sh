#! /bin/bash
#
gfortran -c -Wall components_test.f90
if [ $? -ne 0 ]; then
    echo "Compile error."
  exit
fi
#
gfortran components_test.o $HOME/lib/components.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm components_test.o
#
mv a.out components_test
./components_test > components_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm components_test
#
echo "Normal end of execution."
