#! /bin/bash
#
gfortran -c -Wall -fopenmp fft_openmp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv fft_openmp.o $HOME/lib/fft_openmp.o
#
echo "Normal end of execution."
