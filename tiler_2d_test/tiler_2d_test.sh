#! /bin/bash
#
~/bin/tiler_2d
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
ps2png tiler_2d.ps tiler_2d.png
#
echo "Normal end of execution."
