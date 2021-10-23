### ARRAY DYNAMICS SIMULATION ###
# input:        simulation results
# output:       visualisation of overall activity/methylation range
# description:  fluctuation over time
# author:       HR

library(dplyr)
library(stringr)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

if(!dir.exists(paste0(outfol, 'fluct/'))) {
  dir.create(paste0(outfol, 'fluct/'))
}

print('-------------------------------------')
print('VISUALISATION OF OVERALL FLUCTUATIONS')
print('-------------------------------------')


### INPUT ###
fs = list.files(path = outfol,
                pattern = paste0('_c', c, '_met-', met, '_rep'),
                recursive = F, full.names = T)


### MAIN PART ###
# ----- extract mean activity and methylation -----
ovA = list()
ovM = list()

for (f in 1:length(fs)) {
  
  load(fs[f])
  
  # activity
  times = names(SIMresults$A)
  act = sapply(SIMresults$A, function(x) { length(which(x == 1)) / length(x) })
  names(act) = times
  ovA[[f]] = act
  
  # methylation
  times = names(SIMresults$M)
  meth = sapply(SIMresults$M, function(x) { mean(x) })
  names(meth) = times
  ovM[[f]] = meth
  
  # nm = str_split(fs[f], pattern = '_', simplify = T)[, c(2,3)] %>%
  #   paste(collapse = '_')
  # 
  # names(ovA)[f] = nm
  # names(ovM)[f] = nm
}

cols = heat.colors(length(fs), alpha = .5)

#----- plot activity -----
setwd(paste0(outfol, 'fluct/'))

png(paste0('fluctuations_',lattice,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'_met-',met,'.png'),
    width = 9, height = 9, units = 'in',
    res = 300)

par(mfrow = c(2,1))

plot(ovA[[1]],
     type = 'l',
     xlab = 'time [s]',
     ylab = '<A>',
     ylim = c(0,1),
     lty = 'solid',
     col = cols[1])

if(length(fs) > 1) {
  for (i in 2:length(ovA)) {
    
    lines(ovA[[i]],
          lty = 'solid',
          col = cols[i])
    
  }
}

#----- plot methylation -----

plot(ovM[[1]],
     type = 'l',
     xlab = 'time [s]',
     ylab = '<M>',
     ylim = c(9,30),
     lty = 'solid',
     col = cols[1])

if(length(fs) > 1) {
  for (i in 2:length(ovA)) {
    
    lines(ovM[[i]],
          lty = 'solid',
          col = cols[i])
    
  }
}

dev.off()

### OUTPUT ###

save(ovA, file = paste0('A_c', c, '_met-', met, '.RData'))
save(ovM, file = paste0('M_c', c, '_met-', met, '.RData'))

setwd('~/Documents/SYNMICRO/')

write(paste0('plotting activity flustuations finished sucessfully at ', Sys.time()),
      file = paste0('logs/plotfluctuations_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))


# spectrum(x)
# 
# x = ovA[[1]]
# xPSD = (1/length(x))*abs(fft(x)^2)
# fPSD = seq(0,1.0-1/length(x),by=1/length(x)) / 2*pi
# 
# plot(xPSD~fPSD, type='l')

# ltps = str_extract(names(ovA), pattern = 'J[:digit:].[:digit:]*') %>%
#   str_remove_all(pattern = 'J') %>% 
#   as.numeric()
# LTPS = c('solid', 'dashed', 'dotted')
# names(LTPS) = unique(ltps)
# 
# cols = str_extract(names(ovA), pattern = 'r0[:digit:]*$') %>%
#   str_remove_all(pattern = 'r0') %>% 
#   as.numeric()
# cols[is.na(cols)] = 0.5
# COLS = c('blue', 'red', 'green')
# names(COLS) = unique(cols)

# legend('topright',
#        legend = c(paste0('J = ', unique(ltps)),
#                   paste0('r0 = ', unique(cols))),
#        col = c(rep('black', 3), COLS),
#        lty = c(LTPS, rep('solid', 3)))
# 
# legend('bottomleft',
#        legend = c(paste0('J = ', unique(ltps)),
#                   paste0('r0 = ', unique(cols))),
#        col = c(rep('black', 3), COLS),
#        lty = c(LTPS, rep('solid', 3)))

