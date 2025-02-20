#! /bin/bash
#
gfortran -c -Wall steam_interact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o steam_interact \
  steam_interact.o \
  ~/lib/steam.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm steam_interact.o
#
mv steam_interact ~/bin/steam_interact
#
echo "Normal end of execution."
