### ARRAY DYNAMICS SIMULATION ###
# input:        PSDs for c=0 (simulated and experimental)
# output:       activity fluctuation parameter overview
# description:  visualise the system's activity fluctuation over different parameters
#               and compare with experimental data
# author:       HR


library(dplyr)
library(stringr)
library(lattice)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library(latex2exp)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

spectralRamp = brewer.pal(11, "Spectral") %>%
  colorRampPalette()
spectral5000 = spectralRamp(5000)

theme_set(theme_classic())


### INPUT ###
load('results/PSD_experiments.RData')
# load(snakemake@input[['PSDexp']])

fs = Sys.glob('results/SIMresults/*/PSD/A_c0_met-*.RData')
fs = fs[str_detect(fs, coll(met))]
fs = fs[str_detect(fs, lattice)]
fs = fs[str_detect(fs, "_n")]

print(fs)

### MAIN PART ###
# ----- extract simulation parameters, PSDs and PSD characteristics -----

extractPARAMS = function(f, ret_psd) {
  
  load(f)
  params = str_split(f, pattern = '/', simplify = T)[,3]
  J = str_extract(params, pattern = 'J0.[^_]+') %>%
    str_remove('J')
  
  if (lattice == "Kagome") {
    J = as.numeric(J)
  }
  
  r0 = str_extract(params, pattern = '_r[^_]+') %>%
    str_remove('_r') %>%
    as.numeric()
  
  n = str_extract(params, pattern = 'n[^$]+') %>%
    str_remove('n') %>%
    as.numeric()
  
  psd_params = plyr::ldply(PSDchar$allParams)
  psd_params = apply(psd_params, 2, mean)
  
  if (ret_psd) {
    
    psd = plyr::ldply(PSDchar$PSDs)
    psd = t(psd)
    x = reshape2::melt(psd)
    x$Var2 = NULL
    names(x) = c('freq', 'psd')
    
    x$J = J
    x$r0 = r0
    x$n = n
    x$Smax = psd_params[1]
    x$wmax = psd_params[2]
    x$w12 = psd_params[3]
    
    return(x)
    
  } else {
    
    out = c(J=J, r0=r0, n=n, psd_params)
    return(out)
  }
}


PAR = list()
PSD = list()

for (i in 1:length(fs)) {
  f = fs[i]
  PAR[[i]] = extractPARAMS(f, ret_psd=F)
  PSD[[i]] = extractPARAMS(f, ret_psd=T)
}

PAR = plyr::ldply(PAR)
PSD = plyr::ldply(PSD)

# ----- grid of PSDs -----
no_rs = unique(PSD$r0) %>% length()
no_rs = if(no_rs<4) { 5.5 }

no_Js = unique(PSD$J) %>% length()
no_Js = if(no_Js<4) { 5.5 }

no_ns = unique(PSD$n) %>% length()
no_ns = if(no_ns<4) { 5.5 }

if (lattice == "Square") {
  jj = str_split_fixed(PSD$J, coll("-"), Inf)[,1] %>%
    as.numeric()
  PSD$group = as.factor(PSD$r0*jj)
} else {
  PSD$group = as.factor(PSD$r0*PSD$J)
}

PSD$J = as.factor(PSD$J)

p = ggplot(PSD, aes(x=freq, y=psd, col=group)) +
  geom_smooth(aes(group=group)) +
  scale_color_viridis_d() +
  ylab(TeX("$log_{10}(A_{R}(\\omega))(s)$")) +
  xlab(TeX("$log_{10}(\\frac{\\omega}{2 \\cdot \\pi})(Hz)$")) +
  xlim(c(-3.5, -0.5)) +
  theme(aspect.ratio = 1, legend.position = 'none')

gr = p + 
  #geom_point(aes(x=wmax, y=Smax)) +
  #geom_vline(aes(xintercept=w12), linetype='dashed') +
  facet_grid(J~r0)
cs = p + facet_grid(J~.)
rws = p + facet_grid(.~r0)


pn = ggarrange(ggarrange(rws, NULL, ncol = 2, widths = c(no_rs-1,1)),
          ggarrange(gr, cs, ncol = 2, widths = c(no_rs-3,1)),
          align = 'hv',
          heights = c(1,no_Js-1),
          nrow = 2)


pdf(paste0('results/PSD_simulations_J-r0_', lattice, '_met-', met, '.pdf'),
    height = 16, width = 16)
print(pn)
print(gr)
print(cs)
print(rws)
dev.off()


# ----- heatmaps -----
PAR$group = NULL

# format results
Js = sort(PAR$J) %>% unique()
r0s = sort(PAR$r0) %>% unique()

m = array(NA, dim = c(length(Js), length(r0s), ncol(PAR)-3),
          dimnames = list(Js, r0s, names(PAR)[4:ncol(PAR)]))

print(PAR)

for (c in 1:(ncol(PAR)-3)) {

  for (j.1 in Js) {
    for (r.1 in r0s) {
      
      k = which(PAR$J == j.1 & PAR$r0 == r.1)

      print(dimnames(m)[[1]])
      print(dimnames(m)[[2]])

      if (length(k) > 0) {
        m[which(dimnames(m)[[1]] == j.1),which(dimnames(m)[[2]] == r.1),c] = PAR[k, c+3]
      }
      
    }
  } 
   
}

# plotting

plotMATRIX = function(cnt, param) {
  
  tr = if (met == 'RB+') {
    experimental$allParams$`RB+, QEQE`[names(experimental$allParams$`RB+, QEQE`) == param]
    } else {
      experimental$allParams$`RB-, QEEE, 0um MeAsp`[names(experimental$allParams$`RB-, QEEE, 0um MeAsp`) == param]
    }
  
  h = levelplot(cnt %>% t(),
                col.regions=spectral5000,
                panel = function(...) {
                  panel.fill(col = "white")
                  panel.levelplot(...)
                },
                main=param,
                sub=paste0('experiments: ', tr %>% round(4) %>% as.character()),
                xlab='r0',
                ylab='J')
  
  return(as.grob(h))
}

t = paste0('lattice: ', lattice, ', methylation: ', met)

hm = grid.arrange(plotMATRIX(cnt = m[,,1], param = dimnames(m)[[3]][1]),
             plotMATRIX(cnt = m[,,2], param = dimnames(m)[[3]][2]),
             plotMATRIX(cnt = m[,,3], param = dimnames(m)[[3]][3]),
             nrow=1, top=t)
ggsave(filename = paste0('results/PSD_simulations_J-r0_', lattice, '_met-', met, '_params.pdf'),
       plot = hm, device = cairo_pdf(), height = 4, width = 12,
       dpi = 'retina')


### OUTPUT ###
# write(paste0('aggregation of PSDs finished sucessfully at ', Sys.time()),
#       file = unlist(snakemake@output[['agg_psd']]))

