### ARRAY DYNAMICS SIMULATION ###
# input:        overall acitivity + methylation over time for different concentrations
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

print('----------------------------------------------------------')
print('CALCULATE DOSE-RESPONSE BEHAVIOUR OF ARRAY TO PERTURBATION')
print('----------------------------------------------------------')


### INPUT ###
setwd(paste0(outfol, 'dose-response/'))
getwd()

fs = list.files(path = '../fluct', pattern = 'A_c*',
                recursive = F, full.names = T)
fs.M = list.files(path = '../fluct', pattern = 'M_c*',
                recursive = F, full.names = T)

fs = fs[str_detect(fs, coll(met))]
fs.M = fs.M[str_detect(fs.M, coll(met))]

print(fs)

### MAIN PART ###
# remove zero-concentration
# select correct methylation
# fs = fs[-which(str_detect(fs, 'c0_'))]
fs = fs[which(str_detect(fs, paste0('met-', met)))]

# fs.M = fs.M[-which(str_detect(fs.M, 'c0_'))]
fs.M = fs.M[which(str_detect(fs.M, paste0('met-', met)))]

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


CONS = list(activity = list(),
            methylation = list())

for (f in 1:length(fs)) {
  
  load(fs[f])
  load(fs.M[f])
  
  cntA = lapply(ovA, pick.Steady)
  cntM = lapply(ovM, pick.Steady)
  
  CONS[['activity']][[f]] = cntA %>% unlist() %>% as.numeric()
  CONS[['methylation']][[f]] = cntM %>% unlist() %>% as.numeric()
}


names(CONS[['activity']]) = cs
a = unlist(CONS[['activity']])
a.m = sapply(CONS[['activity']], mean)
a.sd = sapply(CONS[['activity']], sd)
a.sd[is.na(a.sd)] = 0


names(CONS[['methylation']]) = cs
m = unlist(CONS[['methylation']])
m.m = sapply(CONS[['methylation']], mean)
m.sd = sapply(CONS[['methylation']], sd)
m.sd[is.na(m.sd)] = 0

# ----- fit logarithmic model -----
nrep = length(CONS[['activity']][[1]])

dr = try(drm(a~rep(cs, each=nrep), fct = LL.4()))

if(! mode(dr) == 'character') {
  
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
  png(paste0('A_dr_',lattice,'_J',paste(J, collapse = '-'),'_r',r_0,'_X',X,'_c',c,'_met-',met,'.png'),
      width = 9, height = 5, units = 'in',
      res = 300)
  
  plot(log10(cs), a.m,
       ylim = c(0,.6),
       pch = 16,
       ylab = TeX('<A> $\\pm$ SD'), xlab = TeX('$log_{10} \\frac{\\[MeAsp\\]}{K_{off}}$'),
       main = TeX('$f(x) = c + \\frac{d-c}{1 + exp(b \\cdot (log(x) - log(e))}$'),
       sub = paste0('b=', coeff[1], ', c=', coeff[2], ', d=', coeff[3], ', e=', coeff[4]))
  arrows(log10(cs), a.m-a.sd, log10(cs), a.m+a.sd, length=0.05, angle=90, code=3)
  lines(y~log10(x+1e-5), data = fkt.cv)
  
  dev.off()
  
  
}

# ----- fit MWC model -----

fitMWC = function(a, m, cs, k0, k1, Kon, Koff) {
  
  c = rep(cs, each=nrep)
  c = c
  
  f0 = function(Kon, Koff, c) {
    c = c*Koff
    f = k0 - k1*m + log((1 + c/Koff) / (1 + c/Kon))
    return(f)
  }
  
  fit = nls(a~1/(1+exp(N*f0(Kon, Koff, c))),
            start= list(N=10),
            trace = T)
  
  
  # lower = list(N=1, Koff=15, Kon=2800),
  # upper = list(N=500, Koff=25, Kon=3200),
  # algorithm = "port",
  
  print(fit)
  
  N = environment(fit[["m"]][["fitted"]])[["internalPars"]][1]
  # Koff = environment(fit[["m"]][["fitted"]])[["internalPars"]][2]
  # Kon = environment(fit[["m"]][["fitted"]])[["internalPars"]][3]
  
  x = c
  y = 1/(1+exp(N*f0(Kon, Koff, c)))
  
  return(list(N=N, fit=fit, x=x, y=y))
}

mwc.fit = fitMWC(a, m, cs, k0, k1, Kon, Koff)

# plotting 
png(paste0('A_Ns_',lattice,'_J',paste(J, collapse = '-'),'_r',r_0,'_X',X,'_c',c,'_met-',met,'.png'),
    width = 9, height = 5, units = 'in',
    res = 300)

plot(log10(cs), a.m,
     ylim = c(min(a.m, mwc.fit$y)-.1, max(a.m, mwc.fit$y)+.1),
     pch = 16,
     ylab = TeX('<A> $\\pm$ SD'), xlab = TeX('$\\frac{\\[MeAsp\\]}{K_{off}}$'),
     main = TeX('$<A>_{MWC} = \\frac{1}{1 + exp(N \\cdot f_0)}$'),
     sub = paste0('N=', round(mwc.fit$N, 4)))
arrows(log10(cs), a.m-a.sd, log10(cs), a.m+a.sd, length=0.05, angle=90, code=3)

points(log10(mwc.fit$x), mwc.fit$y,
       pch=8, col='firebrick', cex=.5)
k = order(mwc.fit$x)
lines(log10(mwc.fit$x[k]), mwc.fit$y[k], col='firebrick')

legend('topright',
       legend = c('simulated', 'fitted'),
       col = c('black', 'firebrick'), pch = c(16, 8))

dev.off()


### OUTPUT ###

perturb = list(data = data.frame(cs=cs, a.m=a.m, a.sd=a.sd),
               dr = dr,
               MWC = mwc.fit)

save(perturb, file = paste0('A_perturb_met-', met, '.RData'))
setwd('~/Documents/SYNMICRO/')

write(paste0('calculating dose-response curves finished sucessfully at ', Sys.time()),
      file = paste0('logs/calculatedoseresponse_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_X',X,'_c',c,'.txt'))



# nrep = length(CONS[[1]])
# DA = data.frame(L0 = NA, L = NA, da = NA)
# for (r in 1:nrep) {
#   
#   # get all activity changes for all ligand concentration changes
#   # matrix with L0 and L in rows and columns
#   # change from i --> j corresponds to L0 --> L
#   cnt = matrix(ncol = length(cs), nrow = length(cs))
#   colnames(cnt) = cs
#   rownames(cnt) = cs
#   for (i in 1:length(cs)) {
#     for (j in 1:length(cs)) {
#       cnt[i,j] = CONS[[j]][r] - CONS[[i]][r]
#     }
#   }
#   
#   cnt = melt(cnt)
#   names(cnt) = c('L0', 'L', 'da')
#   
#   DA = rbind(DA, cnt)
# }
# 
# DA = na.omit(DA)
# DA = DA[-which(DA$da == 0), ]
# DA = DA[which(DA$L > DA$L0), ]
# 
# MWC = function(L0, L, da, Kon, Koff, a0=0.5, n1=0, n2=100) {
#   
#   Ns = seq(n1, n2, length.out=1000)
#   
#   calcDa = function(N, L0, L, da, Kon, Koff, a0=0.5) {
#     Da = 1 / (1 + (1-a0)*((Koff+L0)^N)*(Kon+L)^N)
#     return(Da)
#   }
#   
#   DAs = rep(NA, length(Ns))
#   for (Ni in 1:length(Ns)) {
#     DAs[Ni] = calcDa(N=Ns[Ni], L0, L, da, Kon, Koff, a0=0.5)
#   }
#   
#   mse = (DAs-da)^2
#   bestN = Ns[which(mse == min(mse))]
#   
#   if (length(bestN) > 1) {
#     bestN = bestN[1]
#   }
#   
#   return(bestN)
# }
# 
# DA$N = NA
# for (r in 1:nrow(DA)) {
#   DA$N[r] = MWC(L0=DA$L0[r], L=DA$L[r], da=DA$da[r],
#                 Kon, Koff)
# }
# 
# density(DA$N) %>% plot
# hist(DA$N, breaks=60)
# 
# # plotting
# png(paste0('A_Ns_met-', met, '.png'),
#     width = 9, height = 5, units = 'in',
#     res = 300)
# hist(DA$N,
#      breaks = 60, freq = F,
#      xlab = 'N', ylab = 'frequency',
#      main = '')
# lines(density(DA$N))
# dev.off()
# 
