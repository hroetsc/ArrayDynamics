### ARRAY DYNAMICS SIMULATION ###
# input:        overall activity over time
# output:       Power Spectral Density
# description:  characterisation of fluctuations using the PSD of <A> over time
# author:       HR

library(dplyr)
library(latex2exp)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

if(!dir.exists(paste0(outfol, 'PSD/'))) {
  dir.create(paste0(outfol, 'PSD/'))
}

### INPUT ###
setwd(paste0(outfol, 'PSD/'))
getwd()

load(paste0('../fluct/A_c0_met-', met, '.RData'))


### MAIN PART ###

# ----- calculate Power Spectral Density
calcPSD = function(meanA) {
  # pick steady state
  meanA.steady = meanA[1000:length(meanA)]
  
  # calculate spectral density
  N = length(meanA.steady)
  x = (1/N)*abs(fft(meanA.steady)^2)
  y = seq(0,1-1/N,1/N)
  
  # smoothen x
  x = smooth(x, kind = '3RSS')
  names(x) = y
  
  return(x)
}

PSDs = lapply(ovA, calcPSD)

# ----- characterise PSD -----

# maximal S_R
# frequency of maximal S_R
# integral of function until S_R decreases

logS = sapply(PSDs, log10)
logw = sapply(PSDs, function(x){
  log10(as.numeric(names(x))/(2*pi))
})

rm = which(!is.finite(logw))
if (length(rm) > 0) {
  logS = logS[-rm, ]
  logw = logw[-rm, ]
}

# only consider 1st half of the spectrum
logS = logS[c(1:ceiling(nrow(logS)/2)),]
logw = logS[c(1:ceiling(nrow(logw)/2)),]

max_S = apply(logS, 2, function(S) {
  S[which(S == max(S))][1]
}) %>%
  unlist()

max_w = logw[apply(logS, 2, function(S){
  which(S == max(S))[1]
}) %>% unlist()]


# integral
tol = 1/abs(max_S) * 0.02

integrals = rep(NA, length(tol))
for (j in 1:length(tol)) {
  i = 1
  while (logS[i+1,j] > (logS[i,j]-tol[j]) & logS[i+1,j] < (logS[i,j]+tol[j])) {
    i = i+1
  }
  integrals[j] = sum(logS[1:i,j])
  names(integrals)[j] = i
}


#----- plotting -----
cols = heat.colors(length(PSDs), alpha = .5)

png(paste0('A_c0_met-', met, '.png'),
    width = 9, height = 5, units = 'in',
    res = 300)

x = PSDs[[1]]

cntx = log10(as.numeric(names(x))/(2*pi))
cnty = log10(x)

ymin = floor(min(cnty))
ymax = ceiling(max(cnty))

rm = which(!is.finite(cntx))
if (length(rm) > 0) {
  cntx = cntx[-rm]
  cnty = cnty[-rm]
}

plot(cntx, cnty,
     type='l',
     ylab = TeX("$log_{10}(A_{R}(\\omega))(s)$"),
     xlab = TeX("$log_{10}(\\frac{\\omega}{2 \\cdot \\pi})(Hz)$"),
     ylim = c(ymin, ymax),
     col = cols[1])

points(max_w[1], max_S[1], pch=16, cex=1.2, col=cols[1])

polygon(x = c(cntx[1],
              cntx[1:as.numeric(names(integrals))[1]],
              cntx[as.numeric(names(integrals))[1]]),
        y = c(ymin,
              cnty[1:as.numeric(names(integrals))[1]],
              ymin),
        col = cols[1], border = NA)

if(length(PSDs) > 1) {
  for (i in 2:length(PSDs)) {
    
    x = PSDs[[i]]
    cntx = log10(as.numeric(names(x))/(2*pi))
    cnty = log10(x)
    
    rm = which(!is.finite(cntx))
    if (length(rm) > 0) {
      cntx = cntx[-rm]
      cnty = cnty[-rm]
    }
    
    lines(cntx, cnty,
         type='l', col = cols[i])
    
    points(max_w[i], max_S[i], pch=16, cex=1.2, col=cols[i])
    
    polygon(x = c(cntx[1],
                  cntx[1:as.numeric(names(integrals))[i]],
                  cntx[as.numeric(names(integrals))[i]]),
            y = c(ymin,
                  cnty[1:as.numeric(names(integrals))[i]],
                  ymin),
            col = cols[i], border = NA)
    
  }
}

dev.off()

### OUTPUT ###

PSDchar = list(PSDs = PSDs)

save(PSDchar, file = paste0('A_c0_met-', met, '.RData'))
setwd('~/Documents/SYNMICRO/')

