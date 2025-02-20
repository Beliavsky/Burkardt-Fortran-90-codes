#! /bin/bash
#
~/bin/frieze > frieze_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
ps2png cells.ps cells.png
ps2png frieze.ps frieze.png
ps2png pattern.ps pattern.png
ps2png region.ps region.png
#
echo "Normal end of execution."
