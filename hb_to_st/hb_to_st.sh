#! /bin/bash
#
gfortran -c -Wall hb_to_st.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran hb_to_st.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hb_to_st.o
#
chmod ugo+x a.out
mv a.out ~/bin/hb_to_st
#
echo "Normal end of execution."
