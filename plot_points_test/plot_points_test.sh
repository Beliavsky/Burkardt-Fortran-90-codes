#! /bin/bash
#
$HOME/bin/plot_points < plot_points_input.txt > plot_points_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."

