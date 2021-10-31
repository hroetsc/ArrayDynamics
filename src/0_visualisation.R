### ARRAY DYNAMICS SIMULATION ###
# input:        simulation results
# output:       GIFs visualising simulated array dynamics
# description:  visualisation of grid activity and methylation states
# author:       HR

library(dplyr)
library(reshape2)
library(lattice)
library(animation)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

if(!dir.exists(paste0(outfol, 'GIFs/'))) {
  dir.create(paste0(outfol, 'GIFs/'))
}

print('--------------------------------')
print('VISUALISATION OF ARRAY BEHAVIOUR')
print('--------------------------------')


### INPUT ###
fs = list.files(path = outfol, full.names = T,
                pattern = paste0('_rep1.RData'))
fs = fs[str_detect(fs, coll(met))]
fs = fs[str_detect(fs, paste0('c', c, '_'))]

print(fs)

load(fs)

setwd(paste0(outfol, 'GIFs/'))

### MAIN PART ###
# ----- convert vector back into matrix -----
Lm = melt(L)
Lm$value = as.character(Lm$value)

transformVector = function(V, Lm) {
  
  Vd = data.frame(value = c(1:n) %>% as.character(),
                  x = V)
  
  mx = left_join(Lm, Vd)
  mx$value = NULL
  MX = dcast(mx, Var1~Var2)[, 2:(X+1)]
  
  return(MX) 
}

# ----- visualise lattice structure -----
# Lg = L
# Lg[!is.na(Lg)] = 1
# lt = levelplot(Lg %>% t(),
#                col.regions = heat.colors(1),
#                xlab = '', ylab = '',
#                colorkey = F,
#                panel = function(...) {
#                  panel.fill(col = "white")
#                  panel.levelplot(...)
#                })
# x11()
# print(lt)
# dev.copy2eps(file = 'results/Kagome_lattice.eps')


# ----- visualise activity/methylation changes -----

K = length(SIMresults$A)

saveGIF(expr = for (k in 1:K) {
  if (k %% 10 == 0) {
    levelplot(transformVector(SIMresults$A[[k]], Lm) %>% t(),
              col.regions = heat.colors(2, rev = T),
              xlab = '', ylab = '',
              colorkey = F,
              scales = list(x = list(at = seq(0,X,10)),
                            y = list(at = seq(0,X,10))),
              panel = function(...) {
                panel.fill(col = "white")
                panel.levelplot(...)
              },
              main = paste0('activity at ', names(SIMresults$A)[k] %>% as.numeric() %>% round(4))) %>%
      print()
  }
},
movie.name = paste0('SIMresultsA_c', c,
                    '_met-', met,
                    '_rep', rep,
                    '.gif'),
img.name = 'SIMresultsA',
convert = 'magick',
clean = T)

gc()


saveGIF(expr = for (k in 1:K) {
  if (k %% 10 == 0) {
    levelplot(transformVector(SIMresults$M[[k]], Lm) %>% t(),
              col.regions = heat.colors(Mtot, rev = T),
              xlab = '', ylab = '',
              colorkey = F,
              scales = list(x = list(at = seq(0,X,10)),
                            y = list(at = seq(0,X,10))),
              panel = function(...) {
                panel.fill(col = "white")
                panel.levelplot(...)
              },
              main = paste0('methylation at ', names(SIMresults$M)[k] %>% as.numeric() %>% round(4))) %>%
      print() 
  }
},
movie.name = paste0('SIMresultsM_c', c,
                    '_met-', met,
                    '_rep', rep,
                    '.gif'),
img.name = 'SIMresultsM',
convert = 'magick',
clean = T)

# ----- increase speed of GIF -----
fps = 10
d = 100/fps

system(paste0('convert -delay ', ceiling(d), ' SIMresultsA_c', c, '_met-', met, '_rep', rep, '.gif', ' SIMresultsA_c', c, '_met-', met, '_rep', rep, '.gif'))
system(paste0('convert -delay ', ceiling(d), ' SIMresultsM_c', c, '_met-', met, '_rep', rep, '.gif', ' SIMresultsM_c', c, '_met-', met, '_rep', rep, '.gif'))

setwd('~/Documents/SYNMICRO/')

### OUTPUT ###
write(paste0('plotting grids finished sucessfully at ', Sys.time()),
      file = paste0('logs/plotgrids_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))
