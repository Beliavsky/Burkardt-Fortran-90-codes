#! /bin/bash
#
~/bin/grf_to_eps heawood.grf heawood.eps
ps2png heawood.eps heawood.png
#
~/bin/grf_to_eps knightstour.grf knightstour.eps
ps2png knightstour.eps knightstour.png
#
~/bin/grf_to_eps petersen.grf petersen.eps
ps2png petersen.eps petersen.png
#
~/bin/grf_to_eps tutte.grf tutte.eps
ps2png tutte.eps tutte.png
#
#
echo "Normal end of execution."
