### ARRAY DYNAMICS SIMULATION ###
# input:        FRET data
# output:       parameters (dose response)
# description:  parse experimental data and extract relevant parameters
# author:       HR

library(dplyr)
library(stringr)
library(latex2exp)

source('src/constants.R')

### INPUT ###
fs = list.files(path = 'data', pattern = '.txt',
                recursive = T, full.names = T)

### MAIN PART ###
# ----- parse data -----
parseData = function(f) {
  tbl = read.table(f, header = F, sep = '\t', skip = 1)
  
  intv = tbl$V1[2] - tbl$V1[1]
  if (any(is.na(tbl$V1))) {
    k = which(is.na(tbl$V1))[1]
    
    tbl$V1[k:nrow(tbl)] = seq(tbl$V1[k-1]+intv, nrow(tbl)*intv, intv)
  }
  
  if (nrow(tbl) %% 2 == 1) {
    tbl$V1 = tbl$V1 - intv
  }
  
  tbl = tbl[, c(1,2)]
  names(tbl) = c('frequency', 'PSD')
  
  # convert FRET ration into activity
  tbl$PSD = tbl$PSD / (lambda^2)
  # mu = tbl$PSD[350:400] %>% mean()  # eps^2/lambda^2
  # tbl$PSD = tbl$PSD - mu+epsilon
  
  # convert to log scale
  tbl$frequency = log10(tbl$frequency / (2*pi))
  tbl$PSD = log10(tbl$PSD)
  
  # cut off for very low and very high frequencies
  m = which((tbl$frequency > -3.5) & (tbl$frequency < -1.5))
  tbl = tbl[m,]
  
  return(tbl)
}


EXP = list()

for (i in 1:length(fs)) {
  f = fs[i]
  
  tbl = parseData(f)
  nm = str_split(f, pattern = coll('/'), simplify = T)
  
  cnt.met = if (str_detect(nm[,3], 'RBplus')) { 'y' } else { 'n' }
  nm = str_split(nm[, 4], '_', simplify = T)
  if (any(nm == 'Mean')) {
    nm = nm[, -which(nm == 'Mean')]
  }
  
  annot = c(receptor = str_remove(nm[2], 's'),
            c = nm[3],
            A = str_extract(nm[4], pattern = 'A-[^t]+'))
  
  EXP[[i]] = list(annot = annot,
                  tbl = tbl)
}

names(EXP) = fs

# ----- characterise PSDs -----

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

# ----- plotting -----

interesting = which(str_detect(names(EXP), 'A-0.5'))
interesting = c(interesting, 8)
names(interesting) = c('RB-, QEEE, 0um MeAsp',
                       'RB-, QEQE, 30um MeAsp',
                       'RB+, QEQE')


png('results/PSD_experiments.png',
    width = 9, height = 5, units = 'in',
    res = 300)

cols = rainbow(length(interesting))
allParams = list()
for (ii in 1:length(interesting)) {
  
  cnt = EXP[[ interesting[ii] ]][['tbl']]
  psdchar = characterisePSD(freq = cnt$frequency, psd = cnt$PSD)
  
  if (ii == 1) {
    plot(PSD~frequency,
         ylab = TeX("$log_{10}(A_{R}(\\omega))(s)$"),
         xlab = TeX("$log_{10}(\\frac{\\omega}{2 \\cdot \\pi})(Hz)$"),
         data = cnt,
         type = 'l', col=cols[ii])
    
  } else {
    lines(PSD~frequency,
          data = cnt,
          col=cols[ii])
  }
  
  points(y=psdchar[1],x=psdchar[2],
         pch = 16,cex=1.2,col=cols[ii])
  abline(v=psdchar[3], lty='dashed', col=cols[ii])
  
  allParams[[ii]] = psdchar
}

legend('topright',
       legend = names(interesting),
       col = cols, cex=.8,
       lty = rep('solid', 3))
dev.off()

names(allParams) = names(interesting)

### OUTPUT ###
# save(allParams, file = 'results/PSD_experiments.RData')
save(allParams, file = unlist(snakemake@output[['PSDexp']]))


