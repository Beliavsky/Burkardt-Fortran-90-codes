#! /bin/bash
#
gfortran -c -Wall snakes_and_ladders_simulation.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv snakes_and_ladders_simulation.o ~/lib/snakes_and_ladders_simulation.o
#
echo "Normal end of execution."
