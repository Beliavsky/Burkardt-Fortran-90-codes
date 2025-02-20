#! /bin/bash
#
$HOME/bin/orbital_decay > orbital_decay.txt
#
gnuplot < orbital_decay_commands.txt
#
echo "Normal end of execution."
