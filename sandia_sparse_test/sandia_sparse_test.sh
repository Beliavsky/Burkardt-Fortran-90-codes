#! /bin/bash
#
gfortran -c -Wall sandia_sparse_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sandia_sparse_test.o $HOME/lib/sandia_sparse.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sandia_sparse_test.o
#
mv a.out sandia_sparse_test
./sandia_sparse_test > sandia_sparse_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sandia_sparse_test
#
echo "Normal end of execution."
