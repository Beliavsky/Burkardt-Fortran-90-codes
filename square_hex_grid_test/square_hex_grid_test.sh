#! /bin/bash
#
gfortran -c -Wall square_hex_grid_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o square_hex_grid_test square_hex_grid_test.o $HOME/lib/square_hex_grid.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm square_hex_grid_test.o
#
./square_hex_grid_test > square_hex_grid_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm square_hex_grid_test
#

gnuplot < hex_grid_1_1_1_commands.txt
gnuplot < hex_grid_2_2_3_commands.txt
gnuplot < hex_grid_3_4_10_commands.txt
gnuplot < hex_grid_4_6_21_commands.txt
gnuplot < hex_grid_5_8_36_commands.txt
gnuplot < hex_grid_6_10_55_commands.txt
gnuplot < hex_grid_7_12_78_commands.txt
gnuplot < hex_grid_8_14_105_commands.txt
gnuplot < hex_grid_9_16_136_commands.txt
gnuplot < hex_grid_10_18_171_commands.txt
gnuplot < hex_grid_11_20_210_commands.txt
gnuplot < hex_grid_21_39_800_commands.txt
#
mv *boundary.txt ../../datasets/square_hex_grid
mv *.png ../../datasets/square_hex_grid
mv *data.txt ../../datasets/square_hex_grid
rm *commands.txt
#
echo "Normal end of execution."
