#! /bin/bash
#
gfortran -c -Wall r83_np.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r83_np.o ~/lib/r83_np.o
#
echo "Normal end of execution."
