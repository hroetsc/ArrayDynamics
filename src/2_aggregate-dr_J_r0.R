### ARRAY DYNAMICS SIMULATION ###
# input:        activity fluctuation upon change in ligand concentration
# output:       dose-response parameter overview
# description:  visualise the system's perturbation response over different parameters
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
fs = Sys.glob('results/SIMresults/*/dose-response/A_perturb_met-*.RData')

fs = fs[str_detect(fs, coll(met))]
fs = fs[str_detect(fs, lattice)]
fs = fs[str_detect(fs, "_n")]

print(fs)


### MAIN PART ###
# ----- extract dose response curves and fitted coefficients -----

extractPARAMS = function(f) {
  
  print(f)
  
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
  
  if (is.data.frame(perturb$data)) {
    tbl = perturb$data
    tbl$J = J
    tbl$r0 = r0
    tbl$n = n
    
    tbl$b = perturb$dr$coefficients[1] %>% as.numeric()
    tbl$c = perturb$dr$coefficients[2] %>% as.numeric()
    tbl$d = perturb$dr$coefficients[3] %>% as.numeric()
    tbl$e = perturb$dr$coefficients[4] %>% as.numeric()
    
    tbl$N = perturb$MWC$N
    
    return(tbl)  
  }
  
}


PAR = list()

for (i in 1:length(fs)) {
  f = fs[i]
  PAR[[i]] = extractPARAMS(f)
}

PAR = plyr::ldply(PAR)


# ----- grid of dose-response curves -----
no_rs = unique(PAR$r0) %>% length()
no_rs = if(no_rs<4) { 5.5 }

no_Js = unique(PAR$J) %>% length()
no_Js = if(no_Js<4) { 5.5 }

no_ns = unique(PAR$n) %>% length()
no_ns = if(no_ns<4) { 5.5 }

if (lattice == "Square") {
  jj = str_split_fixed(PAR$J, coll("-"), Inf)[,1] %>%
    as.numeric()
  PAR$group = as.factor(PAR$r0*jj)
} else {
  PAR$group = as.factor(PAR$r0*PSD$J)
}


PAR$J = as.factor(PAR$J)

p = ggplot(PAR, aes(x=cs, y=a.m, col=group, group=group)) +
  geom_point(aes(x=cs, y=a.m))+
  geom_errorbar(aes(x=cs, ymin=a.m-a.sd, ymax=a.m+a.sd))+
  geom_path(aes(group=group)) +
  scale_x_continuous(trans='log10') +
  scale_color_viridis_d() +
  ylab(TeX('<A> $\\pm$ SD')) +
  xlab(TeX('$\\frac{\\[MeAsp\\]}{K_{off}}$')) +
  theme(aspect.ratio = 1, legend.position = 'none')

gr = p +
  facet_grid(J~r0)
cs = p + facet_grid(J~.)
rws = p + facet_grid(.~r0)


pn = ggarrange(ggarrange(rws, NULL, ncol = 2, widths = c(no_rs-1,1)),
              ggarrange(gr, cs, ncol = 2, widths = c(no_rs-3,1)),
              align = 'hv',
              heights = c(1,no_Js-1),
              nrow = 2)

pdf(paste0('results/dose-response_simulations_J-r0_', lattice, '_met-', met, '.pdf'),
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
ns = sort(PAR$n) %>% unique()

m = array(NA, dim = c(length(Js), length(r0s), ncol(PAR)-6),
          dimnames = list(Js, r0s, names(PAR)[7:ncol(PAR)]))



for (c in 1:(ncol(PAR)-6)) {
  
  for (j.1 in Js) {
    for (r.1 in r0s) {
      
      k = which(PAR$J == j.1 & PAR$r0 == r.1)
      if (length(k) > 0) {
        m[which(dimnames(m)[[1]] == j.1),which(dimnames(m)[[2]] == r.1),c] = PAR[k[1], c+6]
      }
      
    }
  } 
  
}


# plotting
plotMATRIX = function(cnt, param) {
  
  h = levelplot(cnt %>% t(),
                col.regions=spectral5000,
                panel = function(...) {
                  panel.fill(col = "white")
                  panel.levelplot(...)
                },
                main=param,
                xlab='r0',
                ylab='J')
  
  return(as.grob(h))
}

t = paste0('lattice: ', lattice, ', methylation: ', met)

hm = grid.arrange(plotMATRIX(cnt = m[,,1], param = dimnames(m)[[3]][1]),
                  plotMATRIX(cnt = m[,,5], param = dimnames(m)[[3]][5]),
                  nrow=1, top=t)
# hm = plotMATRIX(cnt = m[,,5], param = dimnames(m)[[3]][5])
ggsave(filename = paste0('results/dose-response_simulations_J-r0_', lattice, '_met-', met, '_params.pdf'),
       plot = hm, device = cairo_pdf(), height = 4, width = 12,
       dpi = 'retina')



### OUTPUT ###
# write(paste0('aggregation of dose-response curves finished sucessfully at ', Sys.time()),
#      file = unlist(snakemake@output[['agg_dr']]))
#  
