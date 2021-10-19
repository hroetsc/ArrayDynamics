### ARRAY DYNAMICS SIMULATION ###
# input:         -
# output:       constants
# description:  list all constants used for simulation
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
  make_option(c('-m', '--met'), type = 'character', default = 'y',
              help = 'simulate methylation?'),
  make_option(c('-a', '--kR'), type = 'double', default = kR,
              help = 'methylation rate constant'),
  make_option(c('-b', '--kB'), type = 'double', default = kB,
              help = 'demethylation rate constant'),
  make_option(c('-o', '--Mtot'), type = 'integer', default = Mtot,
              help = 'total number of methylation sites per CSU'),
  make_option(c('-l', '--lattice'), type = 'character', default = 'Kagome',
              help = 'specify lattice structure (script name in src/{l}.R)'),
  make_option(c('-k', '--rep'), type = 'integer', default = 1,
              help = 'simulation replicate ID'),
  make_option(c('-x', '--X'), type = 'integer', default = X,
              help = 'lattice size (1D)')
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

# format coupling energies
if (str_detect(J, '-')) {
  J = str_split(J, pattern = coll('-'), simplify = T) %>%
    as.numeric()
} else {
  J = as.numeric(J)
}

# shorter simulations for dose-response
if (c != 0) {
  TIME = 1.5e03
}

# with/without methylation
if (met == 'n') {
  kR = 0
  kB = 0
  Mtot = 12
}
# express Kon and c in units of Koff
Kon_ = Kon/Koff

# modify lattice size in case of square lattice
if (lattice != 'Kagome') {
  x = X^2-.25*X^2
  X = sqrt(x) %>% ceiling() %>% format(digits = 0) %>% as.numeric()
}

# ----- create output folder -----

outfol = paste0('results/SIMresults/',lattice,'_J',
                paste(J, collapse = '-'),'_r',r_0,'/')
if (!dir.exists(outfol)) {
  dir.create(outfol, recursive = T)
}

save(opt, file = paste0(outfol, 'CMDs.RData'))

# ----- print all current parameters -----

print('-------------------------------------------------')
print('running simulation with the following parameters:')
print('-------------------------------------------------')
print(unlist(opt))
print('-------------------------------------------------')
print('-------------------------------------------------')

