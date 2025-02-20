#! /bin/bash
#
gfortran -c -Wall sandia_cvt_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sandia_cvt_test.o $HOME/lib/sandia_cvt.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sandia_cvt_test.o
#
mv a.out sandia_cvt_test
./sandia_cvt_test > sandia_cvt_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sandia_cvt_test
#
echo "Normal end of execution."
