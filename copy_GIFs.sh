#!/bin/bash

mkdir results/download/fluctuations/
mkdir results/download/dose-response/
mkdir results/download/PSDs/
mkdir results/download/GIFs/

for l in Square; do
	#for J in 0.44 0.46 0.48; do
	for J in 0.42-0.42-0.42 0.44-0.44-0.44 0.46-0.46-0.46 0.48-0.48-0.48; do
		for r in 3 5 10; do
			#for n in 768 972 3072; do
			for n in 1296; do

			cp -rf results/SIMresults/"$l"_J"$J"_r"$r"_n"$n"/GIFs/SIMresultsA_c0_met-RB+_rep10.gif results/download/GIFs/"$l"_J"$J"_r"$r"_n"$n"_A_c0_met-RB+_rep10.gif
			cp -rf results/SIMresults/"$l"_J"$J"_r"$r"_n"$n"/PSD/A_c0_met-RB+.png results/download/PSDs/"$l"_J"$J"_r"$r"_n"$n"_A_c0_met-RB+.png
			cp -rf results/SIMresults/"$l"_J"$J"_r"$r"_n"$n"/PSD/A_c0_met-RB-.png results/download/PSDs/"$l"_J"$J"_r"$r"_n"$n"_A_c0_met-RB-.png

			done
		done
	done
done



cp -rf results/SIMresults/*_n*/fluct/fluctuations_*.png results/download/fluctuations/
cp -rf results/SIMresults/*_n*/dose-response/A_Ns_*.png results/download/dose-response/
cp -rf results/SIMresults/*_n*/dose-response/A_dr_*.png results/download/dose-response/
