#! /bin/bash
#
gfortran -c -Wall maze.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv maze.o ~/lib/maze.o
#
echo "Normal end of execution."
