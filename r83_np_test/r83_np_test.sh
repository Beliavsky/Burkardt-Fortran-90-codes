#! /bin/bash
#
gfortran -c -Wall r83_np_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r83_np_test r83_np_test.o $HOME/lib/r83_np.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r83_np_test.o
#
./r83_np_test > r83_np_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r83_np_test
#
echo "Normal end of execution."
