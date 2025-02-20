#! /bin/bash
#
gfortran -c -Wall normal_dataset.f90
if [ $? -ne 0 ]; then
  echo "Compile error.f90"
  exit
fi
#
gfortran normal_dataset.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm normal_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/normal_dataset
#
echo "Normal end of execution."
