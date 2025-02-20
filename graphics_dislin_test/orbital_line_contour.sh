#! /bin/bash
#
export DISLIN=/usr/local/dislin
export LD_LIBRARY_PATH=$DISLIN:$LD_LIBRARY_PATH
#
gfortran -c -I$DISLIN/gf/real64 -Wall orbital_line_contour.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
$DISLIN/bin/gf95link -r8 orbital_line_contour
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm orbital_line_contour.o
#
./orbital_line_contour > orbital_line_contour.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm orbital_line_contour
#
echo "Normal end of execution."
