#! /bin/bash
#
#  Create a document to simulate interactive user input:
#
cat <<EOF > input.txt
3.0
10
10
10
0.1
crystal_10_10_10.txt
EOF
#
#  Run the program.
#
~/bin/crystal_coordinates < input.txt > crystal_coordinates_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
