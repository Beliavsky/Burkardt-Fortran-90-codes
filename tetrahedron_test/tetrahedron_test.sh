#! /bin/bash
#
gfortran -c -Wall tetrahedron_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran tetrahedron_test.o $HOME/lib/tetrahedron.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm tetrahedron_test.o
#
mv a.out tetrahedron_test
./tetrahedron_test > tetrahedron_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm tetrahedron_test
#
echo "Normal end of execution."
