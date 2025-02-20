#! /bin/bash
#
gfortran -c -Wall reformat.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran reformat.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm reformat.o
#
chmod ugo+x a.out
mv a.out ~/bin/reformat
#
echo "Normal end of execution."
