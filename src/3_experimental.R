### ARRAY DYNAMICS SIMULATION ###
# input:        FRET data
# output:       parameters (dose response)
# description:  parse experimental data and extract relevant parameters
# author:       HR

library(dplyr)
library(stringr)
library(latex2exp)

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
  
  if (nrow(tbl) > 401) {
    tbl = tbl[c(1:401), ]
  }
  
  tbl = tbl[, c(1,2)]
  names(tbl) = c('frequency', 'PSD')
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

characterisePSD = function() {
  
}

# ----- plotting -----

interesting = lapply(EXP, function(x) {
  any(str_detect(x$annot[names(x$annot == 'A')], '0.5'))
}) %>%
  unlist() %>%
  which()


cols = heat.colors(length(interesting))

for (ii in 1:length(interesting)) {
  cnt = EXP[[ interesting[ii] ]][['tbl']]
  
  if (ii == 1) {
    plot(log10(PSD)~log10(frequency/(2*pi)),
         ylab = TeX("$log_{10}(A_{R}(\\omega))(s)$"),
         xlab = TeX("$log_{10}(\\frac{\\omega}{2 \\cdot \\pi})(Hz)$"),
         data = cnt,
         type = 'l', col=cols[ii])
    
  } else {
    lines(log10(PSD)~log10(frequency/(2*pi)),
          data = cnt,
          col=cols[ii])
  }
}


### OUTPUT ###

