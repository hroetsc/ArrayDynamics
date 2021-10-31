### ARRAY DYNAMICS SIMULATION ###
# input:        -
# output:       master table for all simulations
# description:  create overview of all parameter combinations to be tested
# author:       HR


library(dplyr)
library(tidyr)

### INPUT ###

# lattice = 'Kagome'
# met = c('RB+', 'RB-')
# J = c(.4,.42,.44,.46, .48)
# r_0 = c(0.01, 0.1, 0.3, 1, 3, 10)
# c = c(0, 0.1, 0.3, 1, 3, 10)
# rep = 10

lattice = snakemake@params[['lattice']]
met = snakemake@params[['met']]
J = snakemake@params[['J']]
r_0 = snakemake@params[['r_0']]
c = snakemake@params[['c']]
rep = snakemake@params[['rep']]

print('--------------------------------')
print('PARAMETERS')
print(lattice)
print(met)
print(J)
print(r_0)
print(c)
print(rep)
print('--------------------------------')


### MAIN PART ###
# in case lattice is Square --> all combinations of J
if (lattice == 'Square') {
  
  Jexp = expand.grid(J, J, J) %>%
    as.data.frame() %>%
    unique()
  
  J = apply(Jexp, 1, paste, collapse='-')
}

# create all combinations
MASTER = crossing(lattice, met,J,r_0,c,rep) %>%
  as.data.frame()

k = which(MASTER$met == 'RB+' & MASTER$c > 0)
if (length(k) > 0) {
  MASTER = MASTER[-k,]
}

# MASTER$outfol = paste0('results/SIMresults/',MASTER$lattice,'_J',
#                       paste(MASTER$J, sep = '-'),'_r',MASTER$r_0,'/')

print('----------------------------------------------------------')
print(paste0('CREATED MASTER TABLE FOR ', nrow(MASTER)*rep, ' SIMULATIONS'))
print('----------------------------------------------------------')


### OUTPUT ###
# write.csv(MASTER, 'MASTER.csv', row.names = F)
write.csv(MASTER, file = unlist(snakemake@output[['mastertbl']]), row.names = F)

