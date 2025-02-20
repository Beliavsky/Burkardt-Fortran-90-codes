#! /bin/bash
#
gfortran -c -Wall stripack_bench.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o stripack_bench \
  stripack_bench.o \
  ~/lib/stripack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm stripack_bench.o
#
mv stripack_bench ~/bin/stripack_bench
#
echo "Normal end of execution."

