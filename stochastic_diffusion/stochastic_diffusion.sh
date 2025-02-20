#! /bin/bash
#
gfortran -c -Wall stochastic_diffusion.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv stochastic_diffusion.o ~/lib/stochastic_diffusion.o
#
echo "Normal end of execution."
