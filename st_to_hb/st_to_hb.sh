#! /bin/bash
#
gfortran -c -Wall st_to_hb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran st_to_hb.o $HOME/lib/hb_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm st_to_hb.o
#
chmod ugo+x a.out
mv a.out ~/bin/st_to_hb
#
echo "Normal end of execution."
