#! /bin/bash
#
gfortran -c -Wall c_comment.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran c_comment.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm c_comment.o
#
chmod ugo+x a.out
mv a.out ~/bin/c_comment
#
echo "Normal end of execution."
