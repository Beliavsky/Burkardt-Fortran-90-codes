#! /bin/bash
#
gfortran -c -Wall r8ge_np_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ge_np_test r8ge_np_test.o $HOME/lib/r8ge_np.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ge_np_test.o
#
./r8ge_np_test > r8ge_np_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ge_np_test
#
echo "Normal end of execution."
