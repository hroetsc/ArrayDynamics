#!/bin/bash

for J in 0.44 0.46 0.48; do
	for n in 507 768 972 3072 5043 10092; do

		cp -rf results/SIMresults/Kagome_J"$J"_r10_n"$n"/GIFs/SIMresultsA_c0_met-RB+_rep10.gif results/download/GIFs/Kagome_J"$J"_r10_n"$n"_A_c0_met-RB+_rep10.gif
		cp -rf results/SIMresults/Kagome_J"$J"_r10_n"$n"/PSD/A_c0_met-RB+.png results/download/PSDs/Kagome_J"$J"_r10_n"$n"_A_c0_met-RB+.png
		cp -rf results/SIMresults/Kagome_J"$J"_r10_n"$n"/PSD/A_c0_met-RB-.png results/download/PSDs/Kagome_J"$J"_r10_n"$n"_A_c0_met-RB-.png

	done
done
