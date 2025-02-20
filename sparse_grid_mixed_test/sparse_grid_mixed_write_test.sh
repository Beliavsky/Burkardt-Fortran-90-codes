#! /bin/bash
#
gfortran -c -Wall sparse_grid_mixed_write_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sparse_grid_mixed_write_test.o $HOME/lib/sparse_grid_mixed.o $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sparse_grid_mixed_write_test.o
#
mv a.out sparse_grid_mixed_write_test
./sparse_grid_mixed_write_test > sparse_grid_mixed_write_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sparse_grid_mixed_write_test
#
#  Move sparse grid files to dataset directory.
#
mv *_a.txt ../../datasets/sparse_grid_mixed
mv *_b.txt ../../datasets/sparse_grid_mixed
mv *_r.txt ../../datasets/sparse_grid_mixed
mv *_w.txt ../../datasets/sparse_grid_mixed
mv *_x.txt ../../datasets/sparse_grid_mixed
#
echo "Normal end of execution."
