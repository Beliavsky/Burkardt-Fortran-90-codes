#! /bin/bash
#
gfortran -c -Wall pce_burgers.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran pce_burgers.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm pce_burgers.o
#
chmod ugo+x a.out
mv a.out ~/bin/pce_burgers
#
echo "Normal end of execution."
