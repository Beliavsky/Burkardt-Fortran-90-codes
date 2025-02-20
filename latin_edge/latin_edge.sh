#! /bin/bash
#
gfortran -c -Wall latin_edge.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv latin_edge.o ~/lib/latin_edge.o
#
echo "Normal end of execution."
