
import time
from glob import glob
from os.path import join
import random
import string


benchmarks = "benchmarks/"
logs = "logs/"


# ------------------------------ create master table containing all the simulation parameters ------------------------------


rule mastertable:
  params:
    lattice = features['lattice'],
    met = features['met'],
    J = features['J'],
    r_0 = features['r_0'],
    c = features['c'],
    rep = features['replicates']
  output:
    mastertbl = 'MASTER.csv'
  benchmark:
    join(benchmarks, 'mastertable.json')
  log:
    join(logs, 'mastertable.txt')
  script:
    'MASTER.R'


checkpoint check_mastertable:
  input:
    'MASTER.csv'
  output:
    touch('MASTER.done')


# ------------------------------ run the simulations ------------------------------


rule simulation:
  input:
    mastertbl = 'MASTER.csv'
  output:
    sim = 'results/SIMresults/simulation.txt'
  params:
    cores=10,
    rep=features['replicates']
  benchmark:
    join(benchmarks, 'simulation.json')
  log:
    join(logs, 'simulation.txt')
  run:
    m = pd.read_csv('MASTER.csv', sep=',')
    
    for i,j in m.iterrows():
      lattice, met, J, r_0, c = m['lattice'].iloc[i], m['met'].iloc[i], m['J'].iloc[i], m['r_0'].iloc[i], m['c'].iloc[i]
      shell("Rscript src/0_DynamicMC.R --lattice {lattice} --met {met} --J {J} --r0 {r_0} --c {c} --rep {params.rep} --cores {params.cores}")
      
      param = {'lattice': lattice, 'met': met, 'J': J, 'r_0': r_0, 'c': int(c) if c % 1 == 0 else c}
      logfile = "logs/simulation_{lattice}_met-{met}_J{J}_r{r_0}_c{c}.txt".format(**param)
      print('waiting for:', logfile)
      while not os.path.exists(logfile):
        time.sleep(5)
      #shell('until [ -f {logfile} ]; do; sleep 10; done')

    shell('echo "simulation finished sucessfully!" > results/SIMresults/simulation.txt')


checkpoint check_simulation:
  input:
    sim = 'results/SIMresults/simulation.txt'
  output:
    touch('results/SIMresults/simulation.done')


# ------------------------------ visualise simulated array dynamics ------------------------------


rule plotgrids:
  input:
    sim = 'results/SIMresults/simulation.txt',
    mastertbl = 'MASTER.csv'
  output:
    grids = 'results/SIMresults/plotgrids.txt'
  benchmark:
    join(benchmarks, 'plotgrids.json')
  log:
    join(logs, 'plotgrids.txt')
  run:
    m = pd.read_csv('MASTER.csv', sep=',')
    m_red = m[m['c'] == 0]

    for i,j in m_red.iterrows():
      lattice, met, J, r_0, c = m_red['lattice'].iloc[i], m_red['met'].iloc[i], m_red['J'].iloc[i],\
                                m_red['r_0'].iloc[i], m_red['c'].iloc[i]
      shell("Rscript src/0_visualisation.R --lattice {lattice} --met {met} --J {J} --r0 {r_0} --c {c}")

      param = {'lattice': lattice, 'met': met, 'J': J, 'r_0': r_0, 'c': int(c) if c % 1 == 0 else c}
      logfile = 'logs/plotgrids_{lattice}_met-{met}_J{J}_r{r_0}_c{c}.txt'.format(**param)
      print('waiting for:', logfile)
      while not os.path.exists(logfile):
        time.sleep(5)

    shell('echo "plotted all grids sucessfully!" > results/SIMresults/plotgrids.txt')


checkpoint check_plotgrids:
  input:
    grids = 'results/SIMresults/plotgrids.txt'
  output:
    touch('results/SIMresults/plotgrids.done')



rule plotfluctuations:
  input:
    sim = 'results/SIMresults/simulation.txt',
    mastertbl = 'MASTER.csv'
  output:
    fluct = 'results/SIMresults/plotfluctuations.txt'
  benchmark:
    join(benchmarks, 'plotfluctuations.json')
  log:
    join(logs, 'plotfluctuations.txt')
  run:
    m = pd.read_csv('MASTER.csv', sep=',')

    for i,j in m.iterrows():
      lattice, met, J, r_0, c = m['lattice'].iloc[i], m['met'].iloc[i], m['J'].iloc[i], m['r_0'].iloc[i], m['c'].iloc[i]
      shell("Rscript src/0_visualisation2.R --lattice {lattice} --met {met} --J {J} --r0 {r_0} --c {c}")

      param = {'lattice': lattice, 'met': met, 'J': J, 'r_0': r_0, 'c': int(c) if c % 1 == 0 else c}
      logfile = 'logs/plotfluctuations_{lattice}_met-{met}_J{J}_r{r_0}_c{c}.txt'.format(**param)
      print('waiting for:', logfile)
      while not os.path.exists(logfile):
        time.sleep(5)

    shell('echo "plotted all fluctuations sucessfully!" > results/SIMresults/plotfluctuations.txt')


checkpoint check_plotfluctuations:
  input:
    fluct = 'results/SIMresults/plotfluctuations.txt'
  output:
    touch('results/SIMresults/plotfluctuations.done')


# ------------------------------ PSD spectra ------------------------------


rule calculatepsd:
  input:
    fluct = 'results/SIMresults/plotfluctuations.txt',
    mastertbl = 'MASTER.csv'
  output:
    psd = 'results/SIMresults/calculatepsd.txt'
  benchmark:
    join(benchmarks, 'calculatepsd.json')
  log:
    join(logs, 'calculatepsd.txt')
  run:
    m = pd.read_csv('MASTER.csv', sep=',')
    m_red = m[m['c'] == 0]

    for i,j in m_red.iterrows():
      lattice, met, J, r_0, c = m_red['lattice'].iloc[i], m_red['met'].iloc[i], m_red['J'].iloc[i],\
                                m_red['r_0'].iloc[i], m_red['c'].iloc[i]
      shell("Rscript src/1_PSD.R --lattice {lattice} --met {met} --J {J} --r0 {r_0} --c {c}")

      param = {'lattice': lattice, 'met': met, 'J': J, 'r_0': r_0, 'c': int(c) if c % 1 == 0 else c}
      logfile = 'logs/calculatepsd_{lattice}_met-{met}_J{J}_r{r_0}_c{c}.txt'.format(**param)
      print('waiting for:', logfile)
      while not os.path.exists(logfile):
        time.sleep(5)

    shell('echo "calculated all PSDs sucessfully!" > results/SIMresults/calculatepsd.txt')


checkpoint check_calculatepsd:
  input:
    psd = 'results/SIMresults/calculatepsd.txt'
  output:
    touch('results/SIMresults/calculatepsd.done')


rule experimentalpsd:
  output:
    PSDexp = 'results/PSD_experiments.RData'
  benchmark:
    join(benchmarks, 'experimentalpsd.json')
  log:
    join(logs, 'experimentalpsd.txt')
  script:
    '3_experimental.R'


rule aggregatepsd:
  input:
    psd = 'results/SIMresults/calculatepsd.txt',
    PSDexp = 'results/PSD_experiments.RData'
  output:
    agg_psd = 'results/SIMresults/aggregatepsd.txt'
  benchmark:
    join(benchmarks, 'aggregatepsd.json')
  log:
    join(logs, 'aggregatepsd.txt')
  script:
    '2_aggregate-psd.R'

checkpoint check_aggregatepsd:
  input:
    agg_psd = 'results/SIMresults/aggregatepsd.txt'
  output:
    touch('results/SIMresults/aggregatepsd.done')


# ------------------------------ dose-response ------------------------------


rule calculatedoseresponse:
  input:
    fluct = 'results/SIMresults/plotfluctuations.txt',
    mastertbl = 'MASTER.csv'
  output:
    dr = 'results/SIMresults/calculatedoseresponse.txt'
  benchmark:
    join(benchmarks, 'calculatedoseresponse.json')
  log:
    join(logs, 'calculatedoseresponse.txt')
  run:
    m = pd.read_csv('MASTER.csv', sep=',')
    # concentration is set to 0 to avoid redundancy
    # R script lists files with all concentrations

    m_red = m[(m['met'] == 'RB-') & (m['c'] == 0)]
    print(m_red)

    for i,j in m_red.iterrows():
      lattice, met, J, r_0, c = m_red['lattice'].iloc[i], m_red['met'].iloc[i], m_red['J'].iloc[i],\
                                m_red['r_0'].iloc[i], m_red['c'].iloc[i]
      shell("Rscript src/1_dose-response.R --lattice {lattice} --met {met} --J {J} --r0 {r_0} --c {c}")

      param = {'lattice': lattice, 'met': met, 'J': J, 'r_0': r_0, 'c': int(c) if c % 1 == 0 else c}
      logfile = 'logs/calculatedoseresponse_{lattice}_met-{met}_J{J}_r{r_0}_c{c}.txt'.format(**param)
      print('waiting for:', logfile)
      while not os.path.exists(logfile):
        time.sleep(5)

    shell('echo "calculated all dose-response curves sucessfully!" > results/SIMresults/calculatedoseresponse.txt')


checkpoint check_calculatedoseresponse:
  input:
    dr = 'results/SIMresults/calculatedoseresponse.txt'
  output:
    touch('results/SIMresults/calculatedoseresponse.done')


rule aggregatedoseresponse:
  input:
    dr = 'results/SIMresults/calculatedoseresponse.txt'
  output:
    agg_dr = 'results/SIMresults/aggregatedoseresponse.txt'
  benchmark:
    join(benchmarks, 'aggregatedoseresponse.json')
  log:
    join(logs, 'aggregatedoseresponse.txt')
  script:
    '2_aggregate-dr.R'

checkpoint check_aggregatedoseresponse:
  input:
    agg_dr = 'results/SIMresults/aggregatedoseresponse.txt'
  output:
    touch('results/SIMresults/aggregatedoseresponse.done')


# ------------------------------ alternatives ------------------------------

#checkpoint check_mastertable:
#  input:
#    'MASTER.csv'
#  output:
#    touch('MASTER.done')

#class Checkpoint_simulation:
#
#  def __init__(self, pattern):
#    self.pattern = pattern
#
#  def __call__(self, w):
#    global checkpoints
#
#    # wait until master table is fully created
#    checkpoints.check_mastertable.get(**w)
#
#    mastertbl = pd.read_csv('MASTER.csv', sep=',')
#    lattice, met, J, r_0, c = mastertbl['lattice'], mastertbl['met'], mastertbl['J'], mastertbl['r_0'], mastertbl['c']
#
#    pattern = 'results/SIMresults/'+lattice+'_J'+J'_r{r_0}/_c{c}_met-{met}.done'
#
#    return pattern

#checkpoint check_simulation:
#  input:
#    Checkpoint_simulation('results/SIMresults/{lattice}_J{J}_r{r_0}/_c{c}_met-{met}.done')
#  output:
#    touch('results/SIMresults/simulation.done')



