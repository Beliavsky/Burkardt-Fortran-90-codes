#! /bin/bash
#
gfortran -c -Wall ps_write_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o ps_write_test ps_write_test.o $HOME/lib/ps_write.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ps_write_test.o
#
./ps_write_test > ps_write_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ps_write_test
#
ps2png ps_write_test01.ps ps_write_test01.png
ps2png ps_write_test02.ps ps_write_test02.png
ps2png ps_write_test03.ps ps_write_test03.png
ps2png ps_write_test04.ps ps_write_test04.png
ps2png ps_write_test05.ps ps_write_test05.png
ps2png ps_write_test06.ps ps_write_test06.png
ps2png ps_write_test07.ps ps_write_test07.png
ps2png ps_write_test08.ps ps_write_test08.png
ps2png ps_write_test09.eps ps_write_test09.png
ps2png ps_write_test11.eps ps_write_test11.png
ps2png ps_write_test12.eps ps_write_test12.png
#
echo "Normal end of execution."
