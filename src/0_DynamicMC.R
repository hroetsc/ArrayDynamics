### ARRAY DYNAMICS SIMULATION ###
# input:        constants, lattice adjacency list
# output:       simulated receptor dynamics
# description:  simulate switches in activity and methylation events using
#               Dynamic Monte Carlo
# author:       HR

library(dplyr)
library(reshape2)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))

# parameters ro infer: transition rate constant r_0
# coupling energy J
# dose-response curve with different attractant concentrations c

### MAIN PART ###
##### initialise activity states and methylation ######
# could be easier to sample vector directly
# but: to avoid periodicity in random number generator, sample matrix rows
# and flatten to vector

initialiseA = function(X, L) {
  A = matrix(nrow = X, ncol = X)
  
  set.seed(rep)
  for (x in 1:X) {
    A[x, ] = sample(x = c(0,1), size = X, replace = T)
  }
  
  L[which(!is.na(L))] = 1
  A = A*L
  
  A = melt(A) %>%
    na.omit()
  a = A$value
  names(a) = c(1:n)
  
  return(a)
}


initialiseM = function(X, L) {
  M = matrix(nrow = X, ncol = X)
  
  set.seed(rep)
  for (x in 1:X) {
    
    if (met == 'y') {
      M[x, ] = sample(x = seq(0,Mtot), size = X, replace = T)
    } else {
      M[x, ] = rep(Mtot, X)
    }
    
  }
  
  L[which(!is.na(L))] = 1
  M = M*L
  
  M = melt(M) %>%
    na.omit()
  m = M$value
  names(m) = c(1:n)
  
  return(m)
}


##### calculate activity changes #####
calcDeltaA = function(A) {
  DeltaA = -1*(2*A-1)
  return(DeltaA)
}


##### determining coupling energies in a lattice #####

calcJ = function(i, A, adjL, J) {
  
  if (length(adjL) == n) {
    
    k = adjL[[which(names(adjL) == i)]]
    j = 2*J*(2*A[k]-1)
    return(sum(j))
    
  } else {
    
    j = 0
    for (a in 1:length(adjL)) {
      k = adjL[[a]][[which(names(adjL[[a]]) == i)]] %>% na.omit()
      j = j + sum(2*J[a]*(2*A[k]-1))
    }
    return(j)
    
  }
  
} 

##### energy of activity state #####
calcDeltaH = function(r_0, DeltaA, A, M, adjL, k0, k1, c, Kon_, J, idx = c(1:n), ln=NA) {
  
  # free energy - calculate only once in case the concentration is non-zero
  if (is.na(ln)) {
    G = k0 - k1*M[idx] + log((1 + c) / (1 + c/Kon_))
  } else {
    G = k0 - k1*M[idx] + ln
  }
  
  # coupling energy
  # vectorise!
  H_int = rep(NA, length(idx))
  for (i in 1:length(idx)) { H_int[i] = calcJ(idx[i], A, adjL, J) }
  
  # total energy term
  DeltaH = DeltaA[idx]*(G - H_int)
  
  return(DeltaH)
}


########## simulation ##########

SIMULATION = function(adjL, k0, k1, c, Kon_, J, r_0){
  
  set.seed(rep)
  
  A.SIMresults = list()
  M.SIMresults = list()
  
  # ----- initialisation -----
  A = initialiseA(X, L)
  M = initialiseM(X, L)
  DeltaA = calcDeltaA(A)
  
  A.SIMresults[[1]] = A
  names(A.SIMresults)[1] = 0
  
  M.SIMresults[[1]] = M
  names(M.SIMresults)[1] = 0
  
  # ----- calculate initial rates -----
  # activity rates
  ln = log((1 + c) / (1 + c/Kon_))
  DeltaH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c=0, Kon_, J, ln=ln)
  rateA = r_0*exp(-1*DeltaH)
  
  # methylation rates only where methylation/demethylation can happen
  rateM = kR*(1-A)*(M<Mtot) - kB*A*(M>0)
  
  # ----- iteration -----
  pb = txtProgressBar(min = 0, max = TIME, style = 3)
  t = 0
  cnt = 2
  
  while (t < TIME) {
    setTxtProgressBar(pb, t)
    
    # --- sum of all rates ---
    r = c(rateA, abs(rateM))  # a_j in paper
    R = sum(r)  # a_0 in paper
    
    # --- pick time step ---
    u_1 = runif(1, min = epsilon, max = 1-epsilon)
    tau = (1/R) * log(1/u_1)
    
    # --- pick reaction ---
    u_2 = runif(1, min = epsilon, max = 1-epsilon)
    # faster!
    j = which.max(cumsum(r) > u_2*R)
    
    # --- perform reaction and increment time ---
    # decide which reaction (activity or methylation change) happens
    
    if (j <= length(rateA)) {
      
      # change activity for selected site
      A[j] = A[j] + DeltaA[j]
      # update Delta A
      DeltaA[j] = calcDeltaA(A[j])
      
      # --- calculate new rates ---
      # activity rates
      # update the rates of the neighbors
      if (length(adjL) == n) {
        k = adjL[[which(names(adjL) == j)]]
        
        } else {
          
        k = matrix(NA, 3, 2)
        for (a in 1:length(adjL)) {
          k[a, ] = adjL[[a]][[which(names(adjL[[a]]) == j)]]
        }
        k = as.vector(k) %>% na.omit()
        
        }
      
      # spike in ligand at single time point during simulation
      # if (cnt*dt == cT) {
      #   newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c=c, Kon_, J,
      #                     idx = c(j, k))
      # } else {
      #   newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c=0, Kon_, J,
      #                     idx = c(j, k), ln=ln)
      # }
      newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c, Kon_, J,
                        idx = c(j, k), ln=ln)
      rateA[c(j, k)] = r_0*exp(-1*newH)
      
    } else {
      
      j = j-n
      
      # decide whether to methylate or demethylate using the rate's sign
      if (rateM[j] < 0) {
        M[j] = M[j] - 1
      } else if (rateM[j] > 0) {
        M[j] = M[j] + 1
      }
      
      # --- calculate new rates ---
      # activity rates - only of a single position
      # spike in ligand at single time point during simulation
      # if (cnt*dt == cT) {
      #   newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c=c, Kon_, J,
      #                     idx = j)
      # } else {
      #   newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c=0, Kon_, J,
      #                     idx = j, ln=ln)
      # }
      newH = calcDeltaH(r_0, DeltaA, A, M, adjL, k0, k1, c, Kon_, J,
                        idx = j, ln=ln)
      rateA[j] = r_0*exp(-1*newH)
      
    }
    
    # update methylation rates
    rateM[j] = kR*(1-A[j])*(M[j]<Mtot) - kB*A[j]*(M[j]>0)
    
    
    # --- record state of the system every N steps ---
    if (t > cnt*dt) {
      
      A.SIMresults[[cnt]] = A
      names(A.SIMresults)[cnt] = t
      
      M.SIMresults[[cnt]] = M
      names(M.SIMresults)[cnt] = t
      
      cnt = cnt + 1
    }
    
    # --- increment time ---
    t = t + tau
  }
  
  SIMresults = list(A = A.SIMresults,
                    M = M.SIMresults)
  
  return(SIMresults)
}

SIMresults = SIMULATION(adjL, k0, k1, c, Kon_, J, r_0)

### OUTPUT ###
setwd(outfol)
save(SIMresults, file = paste0('_c', c,
                               '_met-', met,
                               '_rep', rep,
                               '.RData'))
setwd('~/Documents/SYNMICRO/')
