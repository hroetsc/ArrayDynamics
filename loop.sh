#!/bin/bash

foo() {

	for m in y n; do
		for c in 0 0.1 0.5 1.0 1.5 5.0 10.0 50.0; do
			for r in 0.1 0.2 0.4 0.6 1.0 1.2 1.5 2.0 5.0 10.0; do
				for J in 0.36 0.38 0.40 0.42 0.44 0.46; do

					Rscript src/0_DynamicMC.R --lattice Kagome --met "$m" --c "$c" --r0 "$r" --J "$J" --rep "$1" ;
					Rscript src/0_visualisation.R --lattice Kagome --met "$m" --c "$c" --r0 "$r" --J "$J" --rep 1;
					Rscript src/0_visualisation2.R --lattice Kagome --met "$m" --c "$c" --r0 "$r" --J "$J" --rep 1;

				done
			done
		done
	done


	for m in y n; do
		for c in 0 0.1 0.5 1.0 1.5 5.0 10.0 50.0; do
			for r in 0.1 0.2 0.4 0.6 1.0 1.2 1.5 2.0 5.0 10.0; do
				for J1 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
					for J2 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
						for J3 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
							Rscript src/0_DynamicMC.R --lattice Square --met "$m" --c "$c" --r0 "$r" --J "$J1"-"$J2"-"$J3" --rep "$1" ;
							Rscript src/0_visualisation.R --lattice Square --met "$m" --c "$c" --r0 "$r" --J "$J1"-"$J2"-"$J3" --rep 1;
							Rscript src/0_visualisation2.R --lattice Square --met "$m" --c "$c" --r0 "$r" --J "$J1"-"$J2"-"$J3" --rep 1;
						done
					done
				done
			done
		done
	done

}


for rep in {1..10}
do
	foo "$rep" &
done


foo2() {

	for m in y n; do
		for r in 0.1 0.2 0.4 0.6 1.0 1.2 1.5 2.0 5.0 10.0; do
			for J in 0.36 0.38 0.40 0.42 0.44 0.46; do
				Rscript src/1_PSD.R --lattice Kagome --met "$m" --c 0 --r0 "$r" --J "$J" --rep 1;
				Rscript src/1_dose-response.R --lattice Kagome --met "$m" --c 0 --r0 "$r" --J "$J" --rep 1;
			done
		done
	done


	for m in y n; do
		for r in 0.1 0.2 0.4 0.6 1.0 1.2 1.5 2.0 5.0 10.0; do
			for J1 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
				for J2 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
					for J3 in 0.36 0.38 0.40 0.42 0.44 0.46 ; do
						Rscript src/1_PSD.R --lattice Square --met "$m" --c 0 --r0 "$r" --J "$J1"-"$J2"-"$J3" --rep 1;
						Rscript src/1_dose-response.R --lattice Square --met "$m" --c 0 --r0 "$r" --J "$J1"-"$J2"-"$J3" --rep 1;
					done
				done
			done
		done
	done

}

foo2