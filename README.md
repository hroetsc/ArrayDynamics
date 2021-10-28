# ArrayDynamics
simulation of E. coli chemoreceptor array dynamics

## General remarks
The pipeline relies on [Conda](https://docs.conda.io/en/latest/) and Snakemake.
In order to install Conda, click on this [link](https://docs.conda.io/en/latest/miniconda.html) and follow the installation guidelines for your respective operating system.  
After installing Conda, you need to install Snakemake. The Snakemake installation procedure is described [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

R, Julia and Python need to be installed. Ideally, install them via Conda. too.

## Simulation parameters
Are specified in `features.yam`. Based on the list, a master table (`MASTER.csv`) is created as first step of the pipeline. The master table contains all parameter combinations tested in the simulation.
