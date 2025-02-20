#! /bin/bash
#
gfortran -c -Wall xerror.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv xerror.o ~/lib/xerror.o
#
echo "Normal end of execution."
