#! /bin/bash
#
gfortran -c -Wall stochastic_diffusion_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o stochastic_diffusion_test \
  stochastic_diffusion_test.o \
  $HOME/lib/stochastic_diffusion.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm stochastic_diffusion_test.o
#
./stochastic_diffusion_test > stochastic_diffusion_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm stochastic_diffusion_test
#
gnuplot < bnt_commands.txt
gnuplot < elman_commands.txt
gnuplot < ntw_commands.txt
gnuplot < xk_commands.txt
#
echo "Normal end of execution."
