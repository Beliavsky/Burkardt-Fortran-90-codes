#! /bin/bash
#
gfortran -c -Wall r8ge_np.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ge_np.o ~/lib/r8ge_np.o
#
echo "Normal end of execution."
