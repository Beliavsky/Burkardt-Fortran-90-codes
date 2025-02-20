#! /bin/bash
#
gfortran -c -Wall maze_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran maze_test.o $HOME/lib/maze.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm maze_test.o
#
mv a.out maze_test
./maze_test > maze_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm maze_test
#
echo "Normal end of execution."
