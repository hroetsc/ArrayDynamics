### ARRAY DYNAMICS SIMULATION ###
# input:        constants, lattice adjacency list
# output:       simulated receptor dynamics
# description:  simulate switches in activity and methylation events using
#               Dynamic Monte Carlo
# author:       HR

library(dplyr)
library(reshape2)
library(foreach)
library(doParallel)

source('src/argparse.R')
source(paste0('src/', lattice, '.R'))


print('------------------------------------------------')
print('DYNAMIC MONTE CARLO SIMULATION OF ARRAY DYNAMICS')
print('------------------------------------------------')


registerDoParallel(cores = cores)

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
  
  for (x in 1:X) {
    
    if (met == 'RB+') {
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
  DeltaA = 1-2*A
  return(DeltaA)
}


##### determining coupling energies in a lattice #####

calcJ = function(i, A, adjL, J) {
  
  if (lattice == 'Kagome') {
    
    k = c(adjL[i,]) %>% na.omit()
    j = 2*J*(2*A[k]-1)
    
    return(sum(j))
    
  } else {
    
    k = adjL[i,]
    
    j = sum(2*J[1]*(2*A[k[1:2] %>% na.omit()]-1)) +
      sum(2*J[2]*(2*A[k[3:4] %>% na.omit()]-1)) +
      sum(2*J[3]*(2*A[k[5:6] %>% na.omit()]-1))
    
    return(j)
    
  }
  
} 

##### energy of activity state #####
calcDeltaH = function(DeltaA, A, M, adjL, k0, k1, c, Kon_, J, idx = c(1:n), ln=NA) {
  
  # free energy - calculate only once in case the concentration is non-zero
  if (is.na(ln)) {
    f = k0 - k1*M[idx] + log((1 + c) / (1 + c/Kon_))
  } else {
    f = k0 - k1*M[idx] + ln
  }
  
  # coupling energy
  # vectorise!
  H_int = rep(NA, length(idx))
  for (i in 1:length(idx)) { H_int[i] = calcJ(idx[i], A, adjL, J) }
  
  # total energy term
  DeltaH = DeltaA[idx]*(f - H_int)
  
  return(DeltaH)
}


########## simulation ##########

SIMULATION = function(adjL, k0, k1, c, Kon_, J, r_0, rep){
  
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
  DeltaH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c, Kon_, J, ln=ln)
  rateA = r_0*exp(-0.5*DeltaH)
  
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
    u_1 = runif(1)
    tau = (1/R) * log(1/u_1)
    
    # --- pick reaction ---
    u_2 = runif(1)
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
      k = adjL[j,] %>%
        na.omit()
      
      # spike in ligand at single time point during simulation
      # if (cnt*dt == cT) {
      #   newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c=c, Kon_, J,
      #                     idx = c(j, k))
      # } else {
      #   newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c=0, Kon_, J,
      #                     idx = c(j, k), ln=ln)
      # }
      newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c, Kon_, J,
                        idx = c(j, k), ln=ln)
      rateA[c(j, k)] = r_0*exp(-0.5*newH)
      
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
      #   newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c=c, Kon_, J,
      #                     idx = j)
      # } else {
      #   newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c=0, Kon_, J,
      #                     idx = j, ln=ln)
      # }
      newH = calcDeltaH(DeltaA, A, M, adjL, k0, k1, c, Kon_, J,
                        idx = j, ln=ln)
      rateA[j] = r_0*exp(-0.5*newH)
      
    }
    
    # update methylation rates
    rateM[j] = kR*(1-A[j])*(M[j]<Mtot) - kB*A[j]*(M[j]>0)
    
    
    # --- record state of the system every N steps ---
    if (t > cnt*dt-dt) {
      
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

tm = proc.time()

foreach(r=1:rep) %dopar% {
  
  SIMresults = SIMULATION(adjL, k0, k1, c, Kon_, J, r_0, r)
  
  ### OUTPUT ###
  save(SIMresults, file = paste0(outfol,
                                 '_c', c,
                                 '_met-', met,
                                 '_rep', r,
                                 '.RData'))
}

elap = proc.time() - tm

stopImplicitCluster()


### OUTPUT ###
write(paste0("R: ",
             round(elap[3], 4), " - r0 ", r_0, ", J ", paste(J, collapse = '-') , ", c ", c, ", met ", met),
      file = paste0(outfol, 'time'), append = T)

write(paste0('simulation finished sucessfully at ', Sys.time()),
      file = paste0('logs/simulation_',lattice,'_met-',met,'_J',paste(J, collapse = '-'),'_r',r_0,'_c',c,'.txt'))

