#! /bin/bash
#
gfortran -c -Wall sandia_sgmgg_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sandia_sgmgg_test.o $HOME/lib/sandia_sgmgg.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sandia_sgmgg_test.o
#
mv a.out sandia_sgmgg_test
./sandia_sgmgg_test > sandia_sgmgg_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sandia_sgmgg_test
#
echo "Normal end of execution."
