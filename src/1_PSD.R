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

print('---------------------------------------------------------')
print('CALCULATE POWER SPECTRAL DENSITY OF ACTIVITY FLUCTUATIONS')
print('---------------------------------------------------------')


### INPUT ###
setwd(paste0(outfol, 'PSD/'))
getwd()

# pick simulations where avg. A is 0.5
# M0=12 RB- <==> c=0
fs = list.files(path = '../fluct', pattern = 'A_c0_',
                recursive = F, full.names = T)
fs = fs[str_detect(fs, coll(met))]
print(fs)

load(fs)


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
  # x = smooth(x, kind = '3RSS')
  
  # convert to log scale
  y = log10(y)
  x = log10(x)
  
  # cut off for very low and very high frequencies
  # m = which((y > -3.5) & (y < -1.5))
  y = y[2:(N/2 +1)]
  x = x[2:(N/2 +1)]
  
  names(x) = y
  return(x)
}

PSDs = lapply(ovA, calcPSD)

# ----- characterise PSD -----

characterisePSD = function(freq, psd) {
  
  # maximum PSD + corresponding frequency
  Smax = max(psd)
  wmax = freq[which(psd == Smax)]
  
  psdnorm = (psd-min(psd)) / (max(psd) - min(psd))
  # psdnorm = (psd-min(psd)) / sum(psd)
  
  # frequency at which half of the integral is reached
  S12 = sum(psdnorm) / 2
  w12 = freq[which(cumsum(psdnorm) >= S12)[1]]
  
  out = c(Smax = Smax,
          wmax = wmax,
          w12 = w12)
  
  return(out)
}

#----- plotting -----
cols = heat.colors(length(PSDs), alpha = .5)

png(paste0('A_c0_met-', met, '.png'),
    width = 9, height = 5, units = 'in',
    res = 300)

allParams = list()
for (ii in 1:length(PSDs)) {
  
  cnt = PSDs[[ii]]
  psdchar = characterisePSD(freq = as.numeric(names(cnt)), psd = cnt)
  
  if (ii == 1) {
    
    plot(y=cnt, x=as.numeric(names(cnt)),
         ylab = TeX("$log_{10}(A_{R}(\\omega))(s)$"),
         xlab = TeX("$log_{10}(\\frac{\\omega}{2 \\cdot \\pi})(Hz)$"),
         ylim = c(min(cnt)-1, max(cnt)+1),
         type = 'l', col=cols[ii])
    
  } else {
    lines(y=cnt, x=as.numeric(names(cnt)),
          col=cols[ii])
  }
  
  points(y=psdchar[1],x=psdchar[2],
         pch = 16,cex=1.2,col=cols[ii])
  abline(v=psdchar[3], lty='dashed', col=cols[ii])
  
  allParams[[ii]] = psdchar
}

dev.off()

### OUTPUT ###

PSDchar = list(PSDs = PSDs,
               allParams = allParams)

save(PSDchar, file = paste0('A_c0_met-', met, '.RData'))
setwd('~/Documents/SYNMICRO/')


write(paste0('calculating PSDs finished sucessfully at ', Sys.time()),
      file = paste0('logs/calculatepsd_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))



