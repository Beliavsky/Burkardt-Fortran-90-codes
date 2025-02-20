#! /bin/bash
#
mpirun -np 3 $HOME/bin/multitask_mpi > multitask_mpi.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."

