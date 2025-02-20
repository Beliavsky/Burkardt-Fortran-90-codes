#! /bin/bash
#
gfortran -c -Wall pdb_to_xyz.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran pdb_to_xyz.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pdb_to_xyz.o
#
chmod ugo+x a.out
mv a.out ~/bin/pdb_to_xyz
#
echo "Normal end of execution."
