#! /bin/bash
#
gfortran -c -Wall matman.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o matman matman.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm matman.o
#
chmod ugo+x matman
mv matman ~/bin/matman
#
echo "Normal end of execution."
