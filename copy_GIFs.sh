#!/bin/bash

for J in 0.4 0.42 0.44 0.46 0.48; do
	for r in 0.01 0.1 0.3 1 3 10; do

		cp -rf results/SIMresults/Kagome_J"$J"_r"$r"/GIFs/SIMresultsA_c0_met-RB+_rep10.gif results/download/GIFs/Kagome_J"$J"_r"$r"_A_c0_met-RB+_rep10.gif

	done
done
