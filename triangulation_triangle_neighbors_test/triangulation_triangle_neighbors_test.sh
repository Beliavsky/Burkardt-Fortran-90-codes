#! /bin/bash
#
if [ -f triangulation_triangle_neighbors_test.txt ]; then
  rm triangulation_triangle_neighbors_test.txt
fi
#
$HOME/bin/triangulation_triangle_neighbors ell > triangulation_triangle_neighbors_test.txt
#
HOME/bin/triangulation_triangle_neighbors lake3 >> triangulation_triangle_neighbors_test.txt
#
echo "Normal end of execution."
