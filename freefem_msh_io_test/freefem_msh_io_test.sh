#! /bin/bash
#
gfortran -c -Wall freefem_msh_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o freefem_msh_io_test freefem_msh_io_test.o $HOME/lib/freefem_msh_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm freefem_msh_io_test.o
#
./freefem_msh_io_test > freefem_msh_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm freefem_msh_io_test
#
echo "Normal end of execution."
