### ARRAY DYNAMICS SIMULATION ###
# input:        constants
# output:       lattice, adjacency matrix
# description:  construct adjacency matrix for Kagome lattice
# author:       HR

library(dplyr)


### MAIN PART ###
##### construct Kagome lattice as matrix #####

getMatrix.Kagome = function(X) {
  L = matrix(nrow = X, ncol = X, byrow = T)
  
  for (x1 in 1:X) {
    
    # for odd row indices, every 4th element is NA
    r = rep(1, X)
    i = 1
    r[1:(i+3) == (i+3)] = NA
    
    # for even row indices, implement shift
    if (x1 %% 2 == 0) {
      r = c(1, NA, r[1:(length(r) - 2)])
    }
    
    L[x1, ] = r
  }
  
  
  
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

L = getMatrix.Kagome(X)

n = which(!is.na(L)) %>% length()
paste0('construct Kagome lattice with ', n, ' elements') %>%
  print()

##### adjacency list of lattice #####
# get all neighbours of a given element

getNeighbors.Kagome = function(L, X) {
  
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
        cntN = rep(NA,4)  # number of nearest neighbors
        
        # odd/even logic is inverted due to zero-padding of L
        # odd column
        if (x1.c %% 2 == 0) {
          
          # first two neighbors are in the same column
          cntN[1] = L[x1.r-1, x1.c]
          cntN[2] = L[x1.r+1, x1.c]
          
          # the other neighbors depend on the position of the NA
          # NA right of current position
          if (!L[x1.r, x1.c-1] %in% c(NA, 0)) {
            cntN[3] = L[x1.r, x1.c-1]
            cntN[4] = L[x1.r+1, x1.c+1]
            
          # NA left of current position
          } else if (!L[x1.r, x1.c+1] %in% c(NA, 0)) {
            cntN[3] = L[x1.r, x1.c+1]
            cntN[4] = L[x1.r+1, x1.c-1]
          }
        
        # even column  
        } else {
          
          cntN[1] = L[x1.r, x1.c+1]
          cntN[2] = L[x1.r, x1.c-1]
          
          cntN[3] = L[x1.r-1, x1.c-1]
          cntN[4] = L[x1.r-1, x1.c+1]
          
        }
        
        cntN = na.omit(cntN)
        
        if (any(cntN == 0)) {
          cntN = cntN[-which(cntN == 0)]
        }
        
        adjL[[i]] = cntN
        names(adjL)[i] = cntE
        i = i + 1
        
      }
      
    }
  }
  
  if (!all(names(adjL) == c(1:n))) {
    print('adjacency list elements are not ordered !!!')
  }
  
  if (! length(adjL) == n) {
    print('!!! adjacency list is incomplete !!!')
  }
  
  return(adjL)
}

adjL = getNeighbors.Kagome(L, X)






