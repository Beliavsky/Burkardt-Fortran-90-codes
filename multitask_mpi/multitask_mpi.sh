#! /bin/bash
#
mpifort -c -Wall multitask_mpi.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mpifort multitask_mpi.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm multitask_mpi.o
mv a.out $HOME/bin/multitask_mpi
#
echo "Normal end of execution."

