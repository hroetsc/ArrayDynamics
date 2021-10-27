### ARRAY DYNAMICS SIMULATION ###
# input:        constants
# output:       lattice, adjacency matrix
# description:  construct adjacency matrix for square lattice
# author:       HR

library(dplyr)


### MAIN PART ###
##### construct square lattice #####

getMatrix.Square = function(X) {
  
  L = matrix(0, nrow = X, ncol = X)
  
  cnt = 1
  for (x1.r in 1:X) {
    for (x1.c in 1:X) {
      
      if (! is.na(L[x1.r, x1.c])) {
        L[x1.r, x1.c] = cnt
        cnt = cnt + 1
      }
      
    }
  }
  
  return(L)
}


L = getMatrix.Square(X)

n = which(!is.na(L)) %>% length()
paste0('construct Square lattice with ', n, ' elements') %>%
  print()

##### adjacency list of lattice #####
# get all neighbours of a given element
# take different kinds of neighbors into account

getNeighbors.Square = function(L, X) {
  
  adjL = list()
  i = 1
  
  # add a frame of 0s to L
  L = cbind(rep(0, X), L, rep(0, X))
  L = rbind(rep(0, X+2), L, rep(0, X+2))
  
  pb = txtProgressBar(min = 2, max = X+1, style = 3)
  
  for (x1.r in 2:(X+1)) {
    setTxtProgressBar(pb, x1.r)
    
    for (x1.c in 2:(X+1)) {
      
      if (! is.na(L[x1.r, x1.c])) {
        
        cntE = L[x1.r, x1.c]
        
        # nearest neighbors
        cntN.1 = rep(NA,2)  # CheW+receptor (vertical)
        cntN.2 = rep(NA,2)  # 2x receptor (horizontal)
        cntN.3 = rep(NA,2)  # 2x CheW (diagonal)
        
        # type 1
        cntN.1[1] = L[x1.r-1, x1.c]
        cntN.1[2] = L[x1.r+1, x1.c]
        
        # type 2
        cntN.2[1] = L[x1.r, x1.c-1]
        cntN.2[2] = L[x1.r, x1.c+1]
        
        # type 3
        cntN.3[1] = L[x1.r-1, x1.c-1]
        cntN.3[2] = L[x1.r+1, x1.c+1]
        
        cntN = c(cntN.1, cntN.2, cntN.3)
        cntN[which(cntN == 0)] = NA
        adjL[[i]] = cntN
        names(adjL)[i] = cntE
        
        i = i + 1 
      }
      
    }
    
  }
  
  
  return(adjL)
}

adj = getNeighbors.Square(L, X)

tmp = plyr::ldply(adj)
rownames(tmp) = tmp$.id
tmp$.id = NULL

tmp = as.matrix(tmp)
adjL = apply(tmp, 2, as.integer) %>%
  as.matrix()

