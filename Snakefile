shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake.utils import min_version

min_version("6.0")
shell.prefix("set -euo pipefail;")

features = yaml.load(open("features.yaml", "r"), Loader=yaml.FullLoader)

snakefiles = "src/"
include: snakefiles + "snakefiles.py"

rule all:
    input:
        mastertbl = 'MASTER.csv',
        sim = 'results/SIMresults/simulation.done',
        grids = 'results/SIMresults/plotgrids.done',
        psd = 'results/SIMresults/calculatepsd.done',
        dr = 'results/SIMresults/calculatedoseresponse.done',
        agg_psd = 'results/SIMresults/aggregatepsd.done',
        agg_dr = 'results/SIMresults/aggregatedoseresponse.done'
        

