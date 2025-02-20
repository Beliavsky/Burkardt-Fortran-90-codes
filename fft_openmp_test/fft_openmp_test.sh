#! /bin/bash
#
gfortran -c -Wall -fopenmp fft_openmp_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -fopenmp -o fft_openmp_test fft_openmp_test.o $HOME/lib/fft_openmp.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm fft_openmp_test.o
#
rm -f fft_openmp_test.txt
#
for threads in 1 2 4 8
do
  echo "Run with "$threads" threads."
  export OMP_NUM_THREADS=$threads
  ./fft_openmp_test >> fft_openmp_test.txt
  if [ $? -ne 0 ]; then
    echo "Run error."
    exit
  fi
done
#
rm fft_openmp_test
#
echo "Normal end of execution."
