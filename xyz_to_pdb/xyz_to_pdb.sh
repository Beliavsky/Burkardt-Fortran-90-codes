#! /bin/bash
#
gfortran -c -Wall xyz_to_pdb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran xyz_to_pdb.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm xyz_to_pdb.o
#
chmod ugo+x a.out
mv a.out ~/bin/xyz_to_pdb
#
echo "Normal end of execution."
