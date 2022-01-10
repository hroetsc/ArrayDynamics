# ArrayDynamics
simulation of E. coli chemoreceptor array dynamics

## General remarks
The pipeline relies on [Conda](https://docs.conda.io/en/latest/) and Snakemake.
In order to install Conda, click on this [link](https://docs.conda.io/en/latest/miniconda.html) and follow the installation guidelines for your respective operating system.  
After installing Conda, you need to install Snakemake. The Snakemake installation procedure is described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

R, Julia and Python need to be installed. Ideally, install them via Conda. too.

Before any execution change the following environment variable: `export JULIA_NUM_THREADS=10`

## Simulation parameters
Are specified in `features.yaml`. Based on the list, a master table (`MASTER.csv`) is created as first step of the pipeline. The master table contains all parameter combinations tested in the simulation.
Make sure to follow the syntax in the template.
The lattice can be either *Kagome* or *Square*. In case of the Kagome lattice, the *J* parameter is a single value; in case of the *Square* lattice, please provide 3 coupling energies for the 3 different link types separated by dashes.

For the RB+ strains, only activities in absence of a chemotactic attractant will be simulated. For the RB- strains, dose-response behaviour to different attractant concentrations is sampled.


