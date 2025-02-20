#! /bin/bash
#
gfortran -c -Wall haar_transform.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv haar_transform.o ~/lib/haar_transform.o
#
echo "Normal end of execution."
