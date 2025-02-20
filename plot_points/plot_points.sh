#! /bin/bash
#
gfortran -c -Wall plot_points.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -g plot_points.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm plot_points.o
#
#  Because you are using a module, you're likely to have a funny
#  compiled module file lying around.
#
rm *.mod
#
chmod ugo+x a.out
mv a.out ~/bin/plot_points
#
echo "Normal end of execution."
