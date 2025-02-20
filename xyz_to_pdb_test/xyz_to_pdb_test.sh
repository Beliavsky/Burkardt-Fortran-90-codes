#! /bin/bash
#
xyz_to_pdb monomer.xyz monomer.pdb > monomer.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
xyz_to_pdb tiny.xyz tiny.pdb > tiny.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
