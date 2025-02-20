#! /bin/bash
#
gfortran -c -Wall hb_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hb_io.o ~/lib/hb_io.o
mv hb_file_module.mod ~/include/hb_file_module.mod
#
echo "Normal end of execution."
