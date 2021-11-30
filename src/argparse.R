### ARRAY DYNAMICS SIMULATION ###
# input:        constants
# output:       simulation parameters
# description:  read and parse command-line arguments + constants
# author:       HR

library(dplyr)
library(optparse)
library(stringr)

source('src/constants.R')

# ----- parse command line arguments -----

option_list = list(
  make_option(c('-j', '--J'), type = 'character', default = J,
              help = 'coupling energy/energies (separate multiple energies by -)'),
  make_option(c('-r', '--r0'), type = 'double', default = r_0,
              help = 'rate constant of activity change'),
  make_option(c('-c', '--c'), type = 'double', default = c,
              help = 'ligand concentration in units of Koff'),
  make_option(c('-m', '--met'), type = 'character', default = 'RB+',
              help = 'simulate methylation?'),
  make_option(c('-a', '--kR'), type = 'double', default = kR,
              help = 'methylation rate constant'),
  make_option(c('-b', '--kB'), type = 'double', default = kB,
              help = 'demethylation rate constant'),
  make_option(c('-o', '--Mtot'), type = 'integer', default = Mtot,
              help = 'total number of methylation sites per CSU'),
  make_option(c('-l', '--lattice'), type = 'character', default = 'Kagome',
              help = 'specify lattice structure (script name in src/{l}.R)'),
  make_option(c('-k', '--rep'), type = 'integer', default = 10,
              help = 'simulation replicate ID'),
  make_option(c('-x', '--X'), type = 'integer', default = X,
              help = 'lattice size (1D)'),
  make_option(c('-n', '--cores'), type = 'integer', default = 10,
              help = 'number of cores for parallel execution of simulation'),
  make_option(c('-t', '--time'), type = 'integer', default = TIME,
              help = 'simulation time (s)')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ----- re-format command line arguments -----
J = opt$J
r_0 = opt$r0
c = opt$c
met = opt$met
kR = opt$kR
kB = opt$kB
Mtot = opt$Mtot
lattice = opt$lattice
rep = opt$rep
X = opt$X
cores = opt$cores
TIME = opt$time

# format coupling energies
if (str_detect(J, '-')) {
  J = str_split(J, pattern = coll('-'), simplify = T) %>%
    as.numeric()
} else {
  J = as.numeric(J)
}
opt$J = J

# shorter simulations for dose-response
if (c != 0) {
  TIME = as.integer(TIME - ceiling(0.25*TIME))
  opt$time = TIME
}

# with/without methylation
if (met == 'RB-') {
  kR = 0
  kB = 0
  Mtot = as.integer(12)

  opt$Mtot = Mtot
  opt$kR = kR
  opt$kB = kB
}

# express Kon and c in units of Koff
Kon_ = Kon/Koff

# modify lattice size in case of square lattice
if (lattice != 'Kagome') {
  # x = X^2-.25*X^2
  # X = sqrt(x) %>% ceiling() %>% format(digits = 0) %>% as.numeric() %>% as.integer()
  # opt$X = X
  
  n = X^2
} else {
  n = 0.75*X^2
}


# ----- create output folder -----

outfol = paste0('results/SIMresults/',lattice,'_J',
                paste(J, collapse = '-'),'_r',r_0,'_n',n,'/')
if (!dir.exists(outfol)) {
  dir.create(outfol, recursive = T)
}

save(opt, file = paste0(outfol, 'CMDs.RData'))

# ----- print all current parameters -----

print('-------------------------------------------------')
print('running with the following parameters:')
print('-------------------------------------------------')
print(unlist(opt))
print('-------------------------------------------------')
print('-------------------------------------------------')

# ----- temporary: compare R and Julia simulations -----

# set.seed(42)
# vu_1 = runif(n=1e07)
# vu_2 = runif(n=1e07)

