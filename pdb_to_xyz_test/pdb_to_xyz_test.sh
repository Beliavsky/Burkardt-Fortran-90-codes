#! /bin/bash
#
pdb_to_xyz monomer.pdb monomer.xyz > monomer.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
pdb_to_xyz tiny.pdb tiny.xyz > tiny.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
