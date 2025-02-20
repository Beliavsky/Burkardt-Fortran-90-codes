#! /bin/bash
#
gfortran -c -Wall sandia_sgmgg.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sandia_sgmgg.o ~/lib/sandia_sgmgg.o
#
echo "Normal end of execution."
