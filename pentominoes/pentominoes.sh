#! /bin/bash
#
mkdir temp
cd temp
rm -f *
~/bin/f90split ../pentominoes.f90
#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -Wall $FILE
  if [ $? -ne 0 ]; then
    echo "Compile error."
    exit
  fi
done
rm *.f90
#
ar qc libpentominoes.a *.o
rm *.o
#
mv libpentominoes.a ~/lib
cd ..
rmdir temp
#
echo "Normal end of execution."
