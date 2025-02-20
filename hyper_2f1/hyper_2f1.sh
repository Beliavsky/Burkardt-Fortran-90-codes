#! /bin/bash
#
gfortran -c -Wall hyper_2f1.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hyper_2f1.o ~/lib/hyper_2f1.o
mv hyper_2f1_module.mod ~/lib/hyper_2f1_module.mod
#
echo "Normal end of execution."
