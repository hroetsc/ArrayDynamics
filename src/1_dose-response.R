### ARRAY DYNAMICS SIMULATION ###
# input:        overall acitivity over time for different concentrations
# output:       dose-response characteristic
# description:  system's response to perturbation (ligand concentration)
#               - logarithmic model
#               - Monod-Wyman-Changeux model (MWC)
# author:       HR

library(dplyr)
library(latex2exp)
library(stringr)
library(reshape2)
library(drc)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

if(!dir.exists(paste0(outfol, 'dose-response/'))) {
  dir.create(paste0(outfol, 'dose-response/'))
}

### INPUT ###
setwd(paste0(outfol, 'dose-response/'))
getwd()

fs = list.files(path = '../fluct', pattern = 'A_c*',
                recursive = F, full.names = T)


### MAIN PART ###
# remove zero-concentration
# select correct methylation
# fs = fs[-which(str_detect(fs, 'c0_'))]
fs = fs[which(str_detect(fs, paste0('met-', met)))]

# ----- get concentrations -----
print(fs)
cs = str_split(fs, '_', simplify = T)[,2] %>%
  str_remove('c') %>% 
  as.numeric()

# ----- load all activities (steady-state) -----
# activity for respective c is mean over steady-state interval
pick.Steady = function(meanA) {
  x = meanA[700:length(meanA)]
  return(mean(x))
}

CONS = list()
for (f in 1:length(fs)) {
  
  load(fs[f])
  
  cntA = lapply(ovA, pick.Steady)
  CONS[[f]] = cntA %>% unlist() %>% as.numeric()
  
}

names(CONS) = cs
a = unlist(CONS)

a.m = sapply(CONS, mean)
a.sd = sapply(CONS, sd)
a.sd[is.na(a.sd)] = 0

# ----- fit logarithmic model -----
dr = drm(a~rep(cs, each=length(CONS[[1]])), fct = LL.4())

fkt = function(coeff, cs) {
  b = coeff[1] %>% as.numeric()
  c = coeff[2] %>% as.numeric()
  d = coeff[3] %>% as.numeric()
  e = coeff[4] %>% as.numeric()
  
  x = seq(min(cs), max(cs), .05)
  y = c + ((d-c)/(1 + exp(b*(log(x) - log(e)))))
  
  return(list(x=x, y=y))
}

coeff = dr$coefficients %>% format(scientific=T) 
fkt.cv = fkt(coeff, cs)

# plotting 
png(paste0('A_dr_met-', met, '.png'),
    width = 9, height = 5, units = 'in',
    res = 300)

plot(cs, a.m,
     ylim = c(0,.6),
     pch = 16,
     ylab = TeX('<A> $\\pm$ SD'), xlab = TeX('$\\frac{\\[MeAsp\\]}{K_{off}}$'),
     main = TeX('$f(x) = c + \\frac{d-c}{1 + exp(b \\cdot (log(x) - log(e))}$'),
     sub = paste0('b=', coeff[1], ', c=', coeff[2], ', d=', coeff[3], ', e=', coeff[4]))
arrows(cs, a.m-a.sd, cs, a.m+a.sd, length=0.05, angle=90, code=3)
lines(y~x, data = fkt.cv)

dev.off()

# ----- fit MWC model -----

nrep = length(CONS[[1]])
DA = data.frame(L0 = NA, L = NA, da = NA)
for (r in 1:nrep) {
  
  # get all activity changes for all ligand concentration changes
  # matrix with L0 and L in rows and columns
  # change from i --> j corresponds to L0 --> L
  cnt = matrix(ncol = length(cs), nrow = length(cs))
  colnames(cnt) = cs
  rownames(cnt) = cs
  for (i in 1:length(cs)) {
    for (j in 1:length(cs)) {
      cnt[i,j] = CONS[[j]][r] - CONS[[i]][r]
    }
  }
  
  cnt = melt(cnt)
  names(cnt) = c('L0', 'L', 'da')
  
  DA = rbind(DA, cnt)
}

DA = na.omit(DA)
DA = DA[-which(DA$da == 0), ]
DA = DA[which(DA$L > DA$L0), ]

MWC = function(L0, L, da, Kon, Koff, a0=0.5, n1=0, n2=100) {
  
  Ns = seq(n1, n2, length.out=1000)
  
  calcDa = function(N, L0, L, da, Kon, Koff, a0=0.5) {
    Da = 1 / (1 + (1-a0)*((Koff+L0)^N)*(Kon+L)^N)
    return(Da)
  }
  
  DAs = rep(NA, length(Ns))
  for (Ni in 1:length(Ns)) {
    DAs[Ni] = calcDa(N=Ns[Ni], L0, L, da, Kon, Koff, a0=0.5)
  }
  
  mse = (DAs-da)^2
  bestN = Ns[which(mse == min(mse))]
  
  if (length(bestN) > 1) {
    bestN = bestN[1]
  }
  
  return(bestN)
}

DA$N = NA
for (r in 1:nrow(DA)) {
  DA$N[r] = MWC(L0=DA$L0[r], L=DA$L[r], da=DA$da[r],
                Kon, Koff)
}

density(DA$N) %>% plot
hist(DA$N, breaks=60)

# plotting
png(paste0('A_Ns_met-', met, '.png'),
    width = 9, height = 5, units = 'in',
    res = 300)
hist(DA$N,
     breaks = 60, freq = F,
     xlab = 'N', ylab = 'frequency',
     main = '')
lines(density(DA$N))
dev.off()

### OUTPUT ###

perturb = list(dr = dr,
               N_MWC = DA)

save(perturb, file = paste0('A_perturb_met-', met, '.RData'))
setwd('~/Documents/SYNMICRO/')
