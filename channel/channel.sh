#! /bin/bash
#
gfortran -c -Wall channel.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran channel.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm channel.o
#
chmod ugo+x a.out
mv a.out ~/bin/channel
#
echo "Normal end of execution."
