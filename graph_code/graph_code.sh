#! /bin/bash
#
gfortran -c -Wall graph_code.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran graph_code.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm graph_code.o
#
chmod ugo+x a.out
mv a.out ~/bin/graph_code
#
echo "Normal end of execution."
