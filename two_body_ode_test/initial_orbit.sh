#! /bin/bash
#
$HOME/bin/initial_orbit > initial_orbit.txt
#
gnuplot < initial_orbit_commands.txt
#
echo "Normal end of execution."
