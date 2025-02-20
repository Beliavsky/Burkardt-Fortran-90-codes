#! /bin/bash
#
$HOME/bin/hypercube_exactness cc_d1_o2 5 > hypercube_exactness_test.txt
$HOME/bin/hypercube_exactness cc_d1_o3 5 >> hypercube_exactness_test.txt
$HOME/bin/hypercube_exactness cc_d2_o3x3 5 >> hypercube_exactness_test.txt
$HOME/bin/hypercube_exactness ccgl_d2_o3x2 5 >> hypercube_exactness_test.txt
#
echo "Normal end of execution."

