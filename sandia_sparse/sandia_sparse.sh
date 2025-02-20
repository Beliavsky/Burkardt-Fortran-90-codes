#! /bin/bash
#
gfortran -c -Wall sandia_sparse.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sandia_sparse.o ~/lib/sandia_sparse.o
#
echo "Normal end of execution."
